#!/usr/bin/env python3
from __future__ import annotations

import streamlit as st

from abi_sauce.services.alignment import build_aligned_columns, pairwise_align
from abi_sauce.services.alignment_waveforms import (
    build_aligned_waveforms,
    raw_channels_and_windows,
)
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.ui.alignment_plot import plot_aligned_traces
from abi_sauce.ui.controls import sample_selector


def _first_trace_asset(sample, assets_view):
    """Pick the first TraceAsset for a Sample (prefers primary)."""
    from abi_sauce.models import TraceAsset

    if sample.primary_asset_id:
        a = assets_view.get(sample.primary_asset_id)
        if isinstance(a, TraceAsset):
            return a
    for aid in sample.asset_ids:
        a = assets_view.get(aid)
        if a.__class__.__name__ == "TraceAsset":
            return a
    return None


def _sample_base_count(sample, assets_view) -> int | None:
    """Robust base count for caption: use the sample’s effective sequence."""
    seq = sample.effective_sequence(assets_view)
    return len(seq) if seq else None


def align_page():
    st.title("Pairwise alignment")
    st.caption(
        "Chromatograms aligned by columns. Choose **Query** and **Reference**, "
        "pick settings, then click **Align** to run."
    )

    sm: SampleManager = st.session_state._samples
    if len(sm.list()) < 2:
        st.info("Upload at least two samples.")
        return

    assets_view = sm.assets_view()  # used both in the form and after submit

    ns = "_align_"
    s = st.session_state
    # Persisted defaults
    s.setdefault(ns + "mode", "globalms")
    s.setdefault(ns + "match", 2.0)
    s.setdefault(ns + "mismatch", -1.0)
    s.setdefault(ns + "gap_open", -10.0)
    s.setdefault(ns + "gap_extend", -0.5)
    s.setdefault(ns + "letter_mode", "axis")  # "axis" (fast) | "overlay"
    s.setdefault(ns + "uniform_x", True)
    s.setdefault(ns + "spp", 21)
    s.setdefault(ns + "force_svg", False)
    s.setdefault(ns + "rangeslider", True)
    s.setdefault(ns + "row_h", 240)
    s.setdefault(ns + "q_id", None)  # Query sample id
    s.setdefault(ns + "r_id", None)  # Reference sample id
    s.setdefault(ns + "has_run", False)

    # ---------- Form: pick everything first ----------
    with st.form("align_form", clear_on_submit=False):
        c1, c2 = st.columns(2)
        with c1:
            # Dropdown (selectbox) for Query
            q_id = sample_selector(sm, label="Query sequence", input_type="selectbox")
        with c2:
            # Dropdown (selectbox) for Reference — no exclusion, can be anything
            r_id = sample_selector(
                sm, label="Reference sequence", input_type="selectbox"
            )

        if q_id and r_id:
            q = sm.get(q_id)
            r = sm.get(r_id)
            q_len = _sample_base_count(q, assets_view)
            r_len = _sample_base_count(r, assets_view)
            q_txt = f"{q.name}: {q_len if q_len is not None else '?'} bp"
            r_txt = f"{r.name}: {r_len if r_len is not None else '?'} bp"
            st.caption(f"**Query** — {q_txt}  •  **Reference** — {r_txt}")

        # Hidden by default; power users can expand
        with st.expander("Alignment & view settings", expanded=False):
            modes = ["globalms", "localms", "globalxx", "localxx"]
            c0, c1, c2, c3, c4 = st.columns([1.2, 1, 1, 1, 1])
            with c0:
                mode_in = st.selectbox("Mode", modes, index=modes.index(s[ns + "mode"]))
            with c1:
                match_in = st.number_input("match", value=float(s[ns + "match"]))
            with c2:
                mismatch_in = st.number_input(
                    "mismatch", value=float(s[ns + "mismatch"])
                )
            with c3:
                gap_open_in = st.number_input(
                    "gap_open", value=float(s[ns + "gap_open"])
                )
            with c4:
                gap_extend_in = st.number_input(
                    "gap_extend", value=float(s[ns + "gap_extend"])
                )

            cL, cU, cS, cR, cH = st.columns([1.6, 1, 1, 1, 1.2])
            with cL:
                letter_choice = st.radio(
                    "Letters",
                    ["axis (fast)", "overlay (per-row)"],
                    index=0 if s[ns + "letter_mode"] == "axis" else 1,
                    horizontal=True,
                )
            with cU:
                uniform_in = st.checkbox(
                    "Uniform x per base", value=bool(s[ns + "uniform_x"])
                )
            with cS:
                spp_in = st.number_input(
                    "Samples/base", 5, 101, value=int(s[ns + "spp"]), step=2
                )
            with cR:
                svg_in = st.checkbox("Force SVG", value=bool(s[ns + "force_svg"]))
            with cH:
                row_h_in = st.slider(
                    "Row height (px)", 180, 480, value=int(s[ns + "row_h"]), step=10
                )

            rs_in = st.toggle("Show rangeslider", value=bool(s[ns + "rangeslider"]))

        submitted = st.form_submit_button(
            "Align", type="primary", use_container_width=True
        )

    # Save selections when submitted
    if submitted:
        s[ns + "q_id"] = q_id
        s[ns + "r_id"] = r_id
        s[ns + "mode"] = mode_in
        s[ns + "match"] = float(match_in)
        s[ns + "mismatch"] = float(mismatch_in)
        s[ns + "gap_open"] = float(gap_open_in)
        s[ns + "gap_extend"] = float(gap_extend_in)
        s[ns + "letter_mode"] = (
            "overlay" if letter_choice.startswith("overlay") else "axis"
        )
        s[ns + "uniform_x"] = bool(uniform_in)
        s[ns + "spp"] = int(spp_in)
        s[ns + "force_svg"] = bool(svg_in)
        s[ns + "rangeslider"] = bool(rs_in)
        s[ns + "row_h"] = int(row_h_in)
        s[ns + "has_run"] = bool(q_id and r_id)

    # Quick utility to swap Q/R and re-run
    c_swap = st.container()
    with c_swap:
        if s.get(ns + "q_id") and s.get(ns + "r_id"):

            def _swap_qr():
                s[ns + "q_id"], s[ns + "r_id"] = s[ns + "r_id"], s[ns + "q_id"]
                s[ns + "has_run"] = True
                st.rerun()

            st.button("Swap Query ↔ Reference", on_click=_swap_qr)

    # If not aligned yet, stop here
    if not s[ns + "has_run"]:
        st.info(
            "Choose **Query** and **Reference**, set options, then click **Align**."
        )
        return

    # ---------- Run alignment with saved selections ----------
    q_id = s[ns + "q_id"]
    r_id = s[ns + "r_id"]
    if not (q_id and r_id):
        st.warning("Select both Query and Reference, then click Align.")
        return

    q = sm.get(q_id)
    r = sm.get(r_id)

    q_trace = _first_trace_asset(q, assets_view)
    r_trace = _first_trace_asset(r, assets_view)
    if not q_trace or not r_trace:
        st.info("Both selected samples need an AB1/trace to show chromatograms.")
        return

    q_seq = (q_trace.sequence or "").upper()
    r_seq = (r_trace.sequence or "").upper()
    if not q_seq or not r_seq:
        st.warning("One or both traces lack a basecalled sequence.")
        return

    with st.spinner("Aligning…"):
        res = pairwise_align(
            q_seq,  # Query
            r_seq,  # Reference
            mode=s[ns + "mode"],
            match=s[ns + "match"],
            mismatch=s[ns + "mismatch"],
            gap_open=s[ns + "gap_open"],
            gap_extend=s[ns + "gap_extend"],
            top_n=1,
        )

    if not res:
        st.error("No alignment produced.")
        return

    r0 = res[0]
    mid = "".join(
        "|" if qa == rb and qa != "-" else " "
        for qa, rb in zip(r0.a_aln, r0.b_aln, strict=False)
    )
    st.markdown(
        f"**Score:** {r0.score:.2f} • **Identity (query vs reference):** {100*r0.identity:.2f}%"
    )
    st.code(f"{r0.a_aln}\n{mid}\n{r0.b_aln}")

    # Channels & windows (per trace)
    q_ch, q_windows, q_len, _, _ = raw_channels_and_windows(q_trace)
    r_ch, r_windows, r_len, _, _ = raw_channels_and_windows(r_trace)
    if q_len == 0 or r_len == 0:
        st.info("Could not compute base windows (missing PLOC/channels).")
        return

    cols = build_aligned_columns(
        r0.a_aln, r0.b_aln
    ).columns  # columns map Query↔Reference indices

    # Cache waveform mapping by trace id + simple column signature + params
    @st.cache_data(show_spinner=False)
    def _cols_sig(pairs: list[tuple[int | None, int | None]]):
        """Cache key: derived signature of alignment columns.
        Keys on (len(pairs), count gaps in Query, count gaps in Reference).
        Invalidation: any change in column structure (new alignment) alters this signature.
        """
        return (
            len(pairs),
            sum(1 for iQ, _ in pairs if iQ is None),
            sum(1 for _, iR in pairs if iR is None),
        )

    @st.cache_data(show_spinner=False)
    def _cached_series(
        trace_id: str,
        sig,
        windows,
        channels,
        for_query: bool,
        uniform: bool,
        spp: int,
    ):
        """Cache key: all arguments (Streamlit caches by value of args).
        Invalidation: changing trace_id, sig (alignment), windows/channels (trace content),
        for_query, uniform, or spp will invalidate and recompute.
        """
        return build_aligned_waveforms(
            columns=cols,
            windows=windows,
            channels=channels,
            for_A=for_query,  # our "A" is Query consistently throughout
            uniform_samples_per_base=uniform,
            samples_per_base=int(spp),
        )

    sig = _cols_sig(cols)
    q_series, _ = _cached_series(
        q_trace.id, sig, q_windows, q_ch, True, s[ns + "uniform_x"], s[ns + "spp"]
    )
    r_series, _ = _cached_series(
        r_trace.id, sig, r_windows, r_ch, False, s[ns + "uniform_x"], s[ns + "spp"]
    )

    fig = plot_aligned_traces(
        columns=cols,
        a_series=q_series,
        a_label=f"Query — {q.name}",
        a_gapped=r0.a_aln,
        a_windows=q_windows,
        a_channels=q_ch,
        b_series=r_series,
        b_label=f"Reference — {r.name}",
        b_gapped=r0.b_aln,
        b_windows=r_windows,
        b_channels=r_ch,
        letter_mode=s[ns + "letter_mode"],
        force_svg=bool(s[ns + "force_svg"]),
        show_rangeslider=bool(s[ns + "rangeslider"]),
        row_height=int(s[ns + "row_h"]),
    )
    st.plotly_chart(fig, use_container_width=True)
