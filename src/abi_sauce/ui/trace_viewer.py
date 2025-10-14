#!/usr/bin/env python3
from __future__ import annotations
from typing import Optional, Dict, Sequence, Tuple
import streamlit as st

from abi_sauce.models import TraceAsset
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.trimming import trim_trace_asset_mott
from abi_sauce.ui.plot_helpers import build_trace_fig
from abi_sauce.ui.defaults import INITIAL_BASE_SPAN, BASE_TICK_EVERY, init_state
from abi_sauce.ui.sequence_viewer import render_sequence_block


def _manager() -> FileManager | None:
    return st.session_state.get("_manager")  # type: ignore[return-value]


def _samples() -> SampleManager | None:
    return st.session_state.get("_samples")  # type: ignore[return-value]


def render_trace_asset(asset: TraceAsset) -> None:
    """
    Render an AB1 trace viewer with trimming & actions + a sequence viewer block.
    """
    st.subheader("Trace (AB1)")

    seq = asset.sequence or ""
    ploc = asset.base_positions or []
    quals = asset.qualities or []
    channels: Dict[str, Sequence[int]] = {
        k: (asset.channels.get(k) or []) for k in ("A", "C", "G", "T")
    }

    if not any(channels.values()):
        st.info("No channel data found in this AB1.")
        # Still show the sequence block if we have one
        if seq:
            st.divider()
            render_sequence_block(
                name=asset.name,
                sequence=seq,
                features=None,
                fasta_text=asset.to_fasta(),
                download_filename=f"{asset.name}.fasta",
                title="Sequence (called)",
                colorize=False,
            )
        return

    # ---- Quick diagnostics / friendly fallbacks ----
    dbg_key = f"abif_dbg_{asset.id}"
    show_dbg = st.sidebar.toggle("Show raw ABIF keys", value=False, key=dbg_key)
    if show_dbg:
        st.sidebar.caption(f"Base order: {asset.meta.get('abif_order', '?')}")
        keys = asset.meta.get("abif_keys") or []
        if keys:
            st.sidebar.code("\n".join(map(str, keys)), language="text")
        else:
            st.sidebar.caption("(no ABIF keys captured)")

    if not ploc:
        st.info("No PLOC2 (peak positions) found — peak hover markers disabled.")
    if not quals:
        st.info("No PHRED qualities detected — quality bars hidden.")

    # ---- State init for chart + trim ----
    prefix = f"traceopts_{asset.id}"
    init_state(prefix)

    show_rangeslider = st.session_state[f"{prefix}_rangeslider"]
    show_grid = st.session_state[f"{prefix}_grid"]
    peak_only_hover = st.session_state[f"{prefix}_hover_peaks_only"]
    err_limit = float(st.session_state[f"{prefix}_err"])
    min_len = int(st.session_state[f"{prefix}_minlen"])
    plot_height = int(st.session_state[f"{prefix}_height"])
    initialized = bool(st.session_state[f"{prefix}_initialized"])

    # ---- Compute Mott trim once (used in chart + sequence highlight) ----
    l_idx = r_idx = None  # 0-based inclusive
    trimmed_seq: Optional[str] = None
    if seq and quals:
        out = trim_trace_asset_mott(seq, quals, error_limit=err_limit, min_len=min_len)
        if out:
            l_idx, r_idx, trimmed_seq = out

    # 1-based inclusive for the figure helper
    seq_len = len(seq) if seq else len(ploc)
    trim_range_1b: Optional[Tuple[int, int]] = None
    if l_idx is not None and r_idx is not None and l_idx <= r_idx:
        trim_range_1b = (int(l_idx) + 1, int(r_idx) + 1)

    # ---- Chart ----
    fig = build_trace_fig(
        channels=channels,
        ploc=ploc,
        seq_len=seq_len,
        quals=quals or None,
        seq=seq or None,
        mode="samples",  # fixed for viewer
        resample_method="max",
        show_trim=True,
        trim_range=trim_range_1b,
        uirevision=f"{asset.id}_chrom_topfirst",
        height=plot_height,
        show_grid=show_grid,
        peak_only_hover=peak_only_hover,
        show_rangeslider=show_rangeslider,
        initial_base_span=INITIAL_BASE_SPAN,
        base_tick_every=BASE_TICK_EVERY,
        x_range=None,
        use_initial_range=not initialized,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.session_state[f"{prefix}_initialized"] = True

    # ---- Options (below chart) ----
    with st.expander("Chart options", expanded=False):
        st.toggle(
            "Show range slider",
            key=f"{prefix}_rangeslider",
            value=st.session_state[f"{prefix}_rangeslider"],
        )
        st.toggle(
            "Show grid", key=f"{prefix}_grid", value=st.session_state[f"{prefix}_grid"]
        )
        st.toggle(
            "Hover only on peaks",
            key=f"{prefix}_hover_peaks_only",
            value=st.session_state[f"{prefix}_hover_peaks_only"],
        )
        st.slider(
            "Plot height (px)",
            320,
            900,
            step=20,
            value=plot_height,
            key=f"{prefix}_height",
            help="Takes effect immediately on toggle/change.",
        )
        c1, c2 = st.columns([1, 1])

        def _reset_view():
            st.session_state[f"{prefix}_initialized"] = False
            st.rerun()

        with c1:
            st.button(
                "Reset chart view", on_click=_reset_view, use_container_width=True
            )
        with c2:
            st.caption("View persists across edits; reset re-applies the default span.")

    # ---- Trim & actions ----
    with st.expander("Trim & actions (Mott)", expanded=False):
        err_key, minlen_key = f"{prefix}_err", f"{prefix}_minlen"
        cur_err = float(st.session_state.get(err_key, err_limit))
        cur_minlen = int(st.session_state.get(minlen_key, min_len))
        with st.form(f"trim_form_{asset.id}", clear_on_submit=False):
            new_err = st.number_input(
                "Error limit (typ. 0.05)",
                0.0,
                0.5,
                value=cur_err,
                step=0.01,
                format="%.3f",
            )
            new_minlen = st.number_input(
                "Minimum kept length", 1, 5000, value=cur_minlen, step=1
            )
            applied = st.form_submit_button("Apply")
        if applied:
            st.session_state[err_key] = max(0.0, min(0.5, float(new_err)))
            st.session_state[minlen_key] = max(1, int(new_minlen))
            st.rerun()

        # Recompute for display below if user changed values inside the expander
        if seq and quals:
            out = trim_trace_asset_mott(
                seq,
                quals,
                error_limit=st.session_state[err_key],
                min_len=st.session_state[minlen_key],
            )
            if out:
                l_idx, r_idx, trimmed_seq = out
                trim_range_1b = (l_idx + 1, r_idx + 1)

        if trimmed_seq:
            fasta = f">{asset.name} | Mott {l_idx+1}-{r_idx+1}\n" + "\n".join(
                trimmed_seq[i : i + 70] for i in range(0, len(trimmed_seq), 70)
            )
            st.download_button(
                "Download Mott-trimmed FASTA",
                data=fasta,
                file_name=f"{asset.name}.mott.fasta",
                use_container_width=True,
            )
            sm = _samples()
            if sm:
                samples_for_asset = [
                    s for s in sm.list() if asset.id in (s.asset_ids or [])
                ]
                if samples_for_asset:
                    st.markdown("##### Apply to Sample")
                    if len(samples_for_asset) == 1:
                        target = samples_for_asset[0]

                        def _apply_one():
                            sm.set_sequence_override(target.id, trimmed_seq)  # type: ignore[arg-type]
                            st.toast(f"Applied trim to “{target.name}”")

                        st.button(f"Apply to “{target.name}”", on_click=_apply_one)
                    else:
                        sel_key = f"trim_target_{asset.id}"
                        names = [s.name for s in samples_for_asset]
                        ids = [s.id for s in samples_for_asset]

                        def _apply_many():
                            sid = st.session_state[sel_key]
                            sm.set_sequence_override(sid, trimmed_seq)  # type: ignore[arg-type]
                            nm = names[ids.index(sid)]
                            st.toast(f"Applied trim to “{nm}”")

                        st.selectbox(
                            "Choose Sample",
                            ids,
                            0,
                            format_func=lambda sid: names[ids.index(sid)],
                            key=sel_key,
                        )
                        st.button("Apply trim", on_click=_apply_many)
                else:
                    st.info(
                        "No Sample references this file yet. Add it to a Sample on the Samples page."
                    )
        else:
            st.caption("(No trimmed sequence available)")

    # ---- Called sequence viewer (with trim highlights) ----
    st.divider()
    if seq:
        # keep_range_0b uses the current (potentially recomputed) indices
        keep_0b: Optional[Tuple[int, int]] = (
            (l_idx, r_idx)
            if (l_idx is not None and r_idx is not None and l_idx <= r_idx)
            else None
        )
        render_sequence_block(
            name=asset.name,
            sequence=seq,
            features=None,
            fasta_text=asset.to_fasta(),
            download_filename=f"{asset.name}.fasta",
            title="Sequence (called)",
            colorize=False,  # plain code block
            plotly_strip=False,  # hide the tiny Plotly strip
        )
    else:
        st.subheader("Sequence (called)")
        st.caption("(No basecalled sequence present in this AB1)")
