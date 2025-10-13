# src/abi_sauce/ui/components.py
from __future__ import annotations
from typing import List, Optional
import streamlit as st

from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.models import Sample

# import trimming helper (used to compute Mott-trim)
from abi_sauce.services.trimming import trim_trace_asset_mott

# import plot helper (single source of truth for plotting)
from abi_sauce.ui.plot_helpers import build_trace_fig

# centralized UI defaults & key init
from abi_sauce.ui.defaults import (
    INITIAL_BASE_SPAN,
    BASE_TICK_EVERY,
    init_state,
)


def asset_table(assets: List[AssetBase]) -> Optional[str]:
    if not assets:
        st.info("No files uploaded yet.")
        return None
    # radio selection
    labels = [f"[{a.kind.value}] {a.name} • {getattr(a, 'length', '')}" for a in assets]
    ids = [a.id for a in assets]
    selected = st.radio("Assets", ids, format_func=lambda x: labels[ids.index(x)])
    return selected


def asset_detail(asset: AssetBase) -> None:
    """
    Renders the detail view for an asset.

    TraceAsset UX:
    - Chart renders FIRST and bigger.
    - All controls (chart options, trim, downloads, sample actions) are BELOW.
    - All bottom expanders start collapsed.
    - X-axis mode is fixed to 'bases' (control hidden).
    """
    # --- Sequence assets ---
    if isinstance(asset, SequenceAsset):
        st.subheader("Sequence")
        st.code(
            asset.sequence[:5000] + ("…" if getattr(asset, "length", 0) > 5000 else ""),
            wrap_lines=True,
        )
        with st.expander("Features", expanded=False):
            if getattr(asset, "features", None):
                st.json(asset.features)
            else:
                st.caption("(No features)")
        fa = asset.to_fasta()
        if fa:
            st.download_button(
                "Download FASTA", data=fa, file_name=f"{asset.name}.fasta"
            )
        return

    # --- Non-trace fallback ---
    if not isinstance(asset, TraceAsset):
        st.write(asset)
        return

    # ========= Trace assets (AB1) =========
    seq = asset.sequence or ""
    ploc = asset.base_positions or []
    quals = asset.qualities or []
    channels = {k: (asset.channels.get(k) or []) for k in ["A", "C", "G", "T"]}

    if not any(channels.values()):
        st.subheader("Trace (AB1)")
        st.info("No channel data found.")
        return

    # ---- Quick diagnostics / friendly fallbacks ----
    # Sidebar micro-toggle shows raw ABIF keys (first 200) and order
    dbg_key = f"abif_dbg_{asset.id}"
    show_dbg = st.sidebar.toggle("Show raw ABIF keys", value=False, key=dbg_key)
    if show_dbg:
        st.sidebar.caption(f"Base order: {asset.meta.get('abif_order', '?')}")
        keys = asset.meta.get("abif_keys") or []
        if keys:
            st.sidebar.code("\n".join(map(str, keys)), language="text")
        else:
            st.sidebar.caption("(no ABIF keys captured)")

    # Inline notices for missing tags
    if not ploc:
        st.info("No PLOC2 (peak positions) found — peak hover markers disabled.")
    if not quals:
        st.info(
            "No PHRED qualities detected (PCON2/letter_annotations) — quality bars hidden."
        )

    # ---- State init for chart + trim ----
    prefix = f"traceopts_{asset.id}"
    init_state(prefix)  # centralizes defaults in one place

    show_rangeslider = st.session_state[f"{prefix}_rangeslider"]
    show_grid = st.session_state[f"{prefix}_grid"]
    peak_only_hover = st.session_state[f"{prefix}_hover_peaks_only"]
    err_limit = st.session_state[f"{prefix}_err"]
    min_len = st.session_state[f"{prefix}_minlen"]
    plot_height = int(st.session_state[f"{prefix}_height"])
    initialized = bool(st.session_state[f"{prefix}_initialized"])

    # ---- Compute trim (used to draw highlight in the chart) ----
    l_idx = r_idx = None
    trimmed_seq = None
    if seq and quals:
        out = trim_trace_asset_mott(seq, quals, error_limit=err_limit, min_len=min_len)
        if out:
            l_idx, r_idx, trimmed_seq = out

    seq_len = len(seq) if seq else len(ploc)
    trim_range = None  # 1-based inclusive
    if l_idx is not None and r_idx is not None and l_idx <= r_idx:
        trim_range = (int(l_idx) + 1, int(r_idx) + 1)

    # ---- ELECTROPHEROGRAM (top) ----
    fig = build_trace_fig(
        channels=channels,
        ploc=ploc,
        seq_len=seq_len,
        quals=quals,
        seq=seq or None,
        mode="samples",
        resample_method="max",
        show_trim=True,
        trim_range=trim_range,
        uirevision=f"{asset.id}_chrom_topfirst",
        height=plot_height,
        show_grid=show_grid,
        peak_only_hover=peak_only_hover,
        show_rangeslider=show_rangeslider,
        initial_base_span=INITIAL_BASE_SPAN,
        base_tick_every=BASE_TICK_EVERY,
        x_range=None,
        use_initial_range=not initialized,  # only first render applies default span
    )
    st.plotly_chart(fig, use_container_width=True)
    st.session_state[f"{prefix}_initialized"] = True

    # ---- Controls & metadata (below chart) ----
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
        colA, colB = st.columns([1, 1])
        with colA:

            def _reset_view():
                st.session_state[f"{prefix}_initialized"] = False
                st.rerun()

            st.button(
                "Reset chart view", on_click=_reset_view, use_container_width=True
            )
        with colB:
            st.caption("View persists across edits; reset re-applies the default span.")

    # --- Trim & actions (debounced form) ---
    with st.expander("Trim & actions (Mott)", expanded=False):
        err_key = f"{prefix}_err"
        minlen_key = f"{prefix}_minlen"
        cur_err = float(st.session_state.get(err_key, 0.05))
        cur_minlen = int(st.session_state.get(minlen_key, 20))
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

        if seq and quals:
            out = trim_trace_asset_mott(
                seq,
                quals,
                error_limit=st.session_state[err_key],
                min_len=st.session_state[minlen_key],
            )
            if out:
                l_idx, r_idx, trimmed_seq = out

        if trimmed_seq:
            fasta = f">{asset.name} | Mott {l_idx+1}-{r_idx+1}\n"
            fasta += "\n".join(
                trimmed_seq[i : i + 70] for i in range(0, len(trimmed_seq), 70)
            )
            st.download_button(
                "Download Mott-trimmed FASTA",
                data=fasta,
                file_name=f"{asset.name}.mott.fasta",
                use_container_width=True,
            )

        sm: Optional[SampleManager] = st.session_state.get("_samples", None)
        if sm and trimmed_seq:
            samples_for_asset = [
                s for s in sm.list() if asset.id in (s.asset_ids or [])
            ]
            if samples_for_asset:
                st.markdown("##### Apply to Sample")
                if len(samples_for_asset) == 1:
                    target = samples_for_asset[0]

                    def _apply_one():
                        sm.set_sequence_override(target.id, trimmed_seq)
                        st.toast(f"Applied trim to “{target.name}”")

                    st.button(f"Apply to “{target.name}”", on_click=_apply_one)
                else:
                    sel_key = f"trim_target_{asset.id}"
                    names = [s.name for s in samples_for_asset]
                    ids = [s.id for s in samples_for_asset]

                    def _apply_many():
                        sid = st.session_state[sel_key]
                        sm.set_sequence_override(sid, trimmed_seq)
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

    with st.expander("Called sequence & downloads", expanded=False):
        if seq:
            st.text(f"Sequence length: {len(seq)} bases")
            fa = asset.to_fasta()
            if fa:
                st.download_button(
                    "Download called sequence (FASTA)",
                    data=fa,
                    file_name=f"{asset.name}.fasta",
                )
        else:
            st.caption("(No basecalled sequence present)")


def sample_table(samples: List[Sample]) -> Optional[str]:
    if not samples:
        st.info("No samples yet — upload files on the Uploads page.")
    labels = [f"{s.name}" for s in samples]
    ids = [s.id for s in samples]
    selected = (
        st.radio("Samples", ids, format_func=lambda x: labels[ids.index(x)])
        if samples
        else None
    )
    return selected


def _rc(seq: str) -> str:
    comp = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return seq.translate(comp)[::-1]


def sample_editor(sample: Sample, fm: FileManager, sm: SampleManager) -> None:
    st.subheader("Sample")

    # --- Name (stateful) ---
    name_key = f"name_input_{sample.id}"
    if name_key not in st.session_state:
        st.session_state[name_key] = sample.name

    colA, colB = st.columns([3, 1])
    with colA:
        st.text_input("Name", key=name_key, label_visibility="visible")
    with colB:

        def _rename():
            sm.rename(sample.id, st.session_state[name_key])
            st.toast("Renamed")

        st.button("Rename", on_click=_rename)

    # --- Primary asset selector (stateful + callback) ---
    if len(sample.asset_ids) > 1:
        st.caption("Primary asset (used for default sequence/features):")
        radio_key = f"primary_{sample.id}"
        if radio_key not in st.session_state:
            if sample.primary_asset_id in sample.asset_ids:
                st.session_state[radio_key] = sample.primary_asset_id
            else:
                st.session_state[radio_key] = sample.asset_ids[0]

        def _set_primary():
            sm.set_primary(sample.id, st.session_state[radio_key])

        st.radio(
            "Primary asset",
            sample.asset_ids,
            index=sample.asset_ids.index(st.session_state[radio_key]),
            label_visibility="collapsed",
            format_func=lambda aid: f"{fm.get(aid).kind.value} • {fm.get(aid).name}",
            key=radio_key,
            on_change=_set_primary,
        )

    # --- Sequence editor (stateful) ---
    seq_default = sample.sequence_override or (
        sample.effective_sequence(fm._assets) or ""
    )
    seq_key = f"seq_buf_{sample.id}"
    if seq_key not in st.session_state:
        st.session_state[seq_key] = seq_default

    st.subheader("Sequence (editable)")
    st.caption("Edits are stored on the sample; original file remains unchanged.")
    st.text_area("sequence", key=seq_key, height=150, label_visibility="collapsed")

    def _save_seq():
        sm.set_sequence_override(sample.id, st.session_state[seq_key])
        st.toast("Saved sequence")

    def _rc_seq():
        s = st.session_state[seq_key]
        comp = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
        rc = s.translate(comp)[::-1]
        st.session_state[seq_key] = rc
        sm.set_sequence_override(sample.id, rc)
        st.toast("Reverse-complemented")

    c3, c4, c5 = st.columns(3)
    with c3:
        st.button("Save sequence", on_click=_save_seq)
    with c4:
        st.button("Reverse complement", on_click=_rc_seq)
    with c5:
        fa = sm.fasta(sample.id)
        if fa:
            st.download_button(
                "Download FASTA", data=fa, file_name=f"{sample.name}.fasta"
            )

    # --- Features editor (stateful) ---
    st.subheader("Features (editable)")
    feats = (
        sample.feature_overrides
        if sample.feature_overrides is not None
        else sample.effective_features(fm._assets)
    )
    edited = st.data_editor(
        feats, num_rows="dynamic", width="stretch", key=f"fe_{sample.id}"
    )
    if st.button("Save features"):
        sm.set_feature_overrides(sample.id, edited)
        st.toast("Saved features")
