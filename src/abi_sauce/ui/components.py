# src/abi_sauce/ui/components.py
from __future__ import annotations
from typing import List, Optional, Dict
import streamlit as st

from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset, AssetKind
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.models import Sample

# import trimming helper (used to compute Mott-trim)
from abi_sauce.services.trimming import trim_trace_asset_mott

# import plot helper (new single-source-of-truth for plotting)
from abi_sauce.ui.plot_helpers import build_trace_fig


def asset_table(assets: List[AssetBase]) -> Optional[str]:
    if not assets:
        st.info("No files uploaded yet.")
        return None
    # radio selection
    labels = [f"[a.kind] {a.name} • {getattr(a, 'length', '')}" for a in assets]
    ids = [a.id for a in assets]
    selected = st.radio("Assets", ids, format_func=lambda x: labels[ids.index(x)])
    return selected


def asset_detail(asset: AssetBase) -> None:
    """
    Renders the detail view for an asset. For TraceAsset this shows the chromatogram,
    and it delegates plot construction to abi_sauce.ui.plot_helpers.build_trace_fig.
    """
    if isinstance(asset, SequenceAsset):
        st.subheader("Sequence")
        st.code(asset.sequence[:5000] + ("…" if asset.length > 5000 else ""))
        with st.expander("Features"):
            st.json(asset.features or {}) if asset.features else st.caption("(No features)")
        fa = asset.to_fasta()
        st.download_button("Download FASTA", data=fa, file_name=f"{asset.name}.fasta")
        return

    if not isinstance(asset, TraceAsset):
        st.write(asset)
        return

    # ---- Trace UI ----
    st.subheader("Trace (AB1)")
    if asset.sequence:
        st.text(f"Sequence length: {len(asset.sequence)} bases")
        fa = asset.to_fasta()
        if fa:
            st.download_button("Download called sequence (FASTA)", data=fa, file_name=f"{asset.name}.fasta")

    if not any(asset.channels.values()):
        st.info("No channel data found.")
        return

    # Prepare channel arrays (ensure lists, default empty)
    channels: Dict[str, List[float]] = {
        "A": asset.channels.get("A") or [],
        "C": asset.channels.get("C") or [],
        "G": asset.channels.get("G") or [],
        "T": asset.channels.get("T") or [],
    }
    ploc = asset.base_positions or []
    quals = asset.qualities or []

    # trimming controls
    with st.expander("Trimming (Mott algorithm)", expanded=True):
        c1, c2, c3 = st.columns([1, 1, 2])
        with c1:
            err_lim = st.slider("Error limit (prob.)", 0.0, 0.20, 0.05, 0.005)
        with c2:
            min_len = st.number_input("Min length", value=20, min_value=1, step=1)
        with c3:
            show_trim = st.checkbox("Show trimmed region", value=True)

    # run trimming (returns 0-based inclusive indexes and trimmed sequence)
    trim = trim_trace_asset_mott(asset.sequence, quals, error_limit=err_lim, min_len=min_len)
    l_idx = r_idx = None
    trimmed_seq = None
    if trim:
        l_idx, r_idx, trimmed_seq = trim

    # Plot mode controls
    c_mode, c_resample = st.columns([1, 1])
    with c_mode:
        mode = st.selectbox("X-axis mode", ["samples", "bases"], help="samples = raw signal sample index with base-number ticks; bases = resampled signals at base positions", index=0)
    with c_resample:
        resample_method = st.selectbox("Resample method (when plotting per-base)", ["max", "mean", "median"], index=0)

    # Build the figure using the new helper
    seq_len = len(asset.sequence) if asset.sequence else len(ploc)
    # Convert trimming indexes (0-based inclusive) to helper's 1-based inclusive trim_range if present
    trim_range = None
    if l_idx is not None and r_idx is not None and l_idx <= r_idx:
        trim_range = (int(l_idx) + 1, int(r_idx) + 1)

    fig = build_trace_fig(
        channels=channels,
        ploc=ploc,
        seq_len=seq_len,
        quals=quals,
        seq=asset.sequence,
        mode=mode,
        resample_method=resample_method,
        show_trim=show_trim,
        trim_range=trim_range,
        uirevision=f"{asset.id}_chrom_v2",
        height=360,
        show_grid=True,
    )

    # Use Streamlit's new width parameter instead of deprecated use_container_width
    st.plotly_chart(fig, width="stretch")

    # Download Mott-trimmed FASTA if trimming produced a result
    if trimmed_seq:
        fasta = f">{asset.name}_mott_trim\n" + "\n".join(
            trimmed_seq[i:i+70] for i in range(0, len(trimmed_seq), 70)
        )
        st.download_button("Download Mott-trimmed FASTA", data=fasta, file_name=f"{asset.name}.mott.fasta")


def sample_table(samples: List[Sample]) -> Optional[str]:
    if not samples:
        st.info("No samples yet — upload files on the Uploads page.")
        return None
    labels = [f"{s.name}" for s in samples]
    ids = [s.id for s in samples]
    selected = st.radio("Samples", ids, format_func=lambda x: labels[ids.index(x)])
    return selected


def _rc(seq: str) -> str:
    comp = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return seq.translate(comp)[::-1]


def sample_editor(sample: Sample, fm: FileManager, sm: SampleManager) -> None:
    st.subheader("Sample")
    c1, c2 = st.columns([3, 1])
    with c1:
        new_name = st.text_input("Name", value=sample.name, label_visibility="visible")
    with c2:
        if st.button("Rename"):
            sm.rename(sample.id, new_name)
            st.toast("Renamed")

    # primary asset selector
    if len(sample.asset_ids) > 1:
        st.caption("Primary asset (used for default sequence/features):")
        st.radio(
            "Primary asset",
            sample.asset_ids,
            index=max(0, sample.asset_ids.index(sample.primary_asset_id) if sample.primary_asset_id in sample.asset_ids else 0),
            label_visibility="collapsed",
            format_func=lambda aid: f"{fm.get(aid).kind.value} • {fm.get(aid).name}",
            key=f"primary_{sample.id}",
            on_change=lambda: sm.set_primary(sample.id, st.session_state[f'primary_{sample.id}'])
        )

    assets_map: Dict[str, AssetBase] = {aid: fm.get(aid) for aid in sample.asset_ids}
    seq = sample.effective_sequence(fm._assets) or ""
    st.subheader("Sequence (editable)")
    st.caption("Edits are stored on the sample; original file remains unchanged.")
    s = st.text_area("sequence", value=(sample.sequence_override or seq), height=150, label_visibility="collapsed")
    c3, c4, c5 = st.columns(3)
    with c3:
        if st.button("Save sequence"):
            sm.set_sequence_override(sample.id, s)
            st.toast("Saved sequence")
    with c4:
        if st.button("Reverse complement"):
            st.session_state[f"seq_buf_{sample.id}"] = _rc(s)
            sm.set_sequence_override(sample.id, _rc(s))
            st.rerun()
    with c5:
        fa = sm.fasta(sample.id)
        if fa:
            st.download_button("Download FASTA", data=fa, file_name=f"{sample.name}.fasta")

    # features (simple table)
    st.subheader("Features (editable)")
    feats = sample.feature_overrides if sample.feature_overrides is not None else sample.effective_features(fm._assets)
    edited = st.data_editor(
        feats,
        num_rows="dynamic",
        width="stretch",
        key=f"fe_{sample.id}"
    )
    if st.button("Save features"):
        sm.set_feature_overrides(sample.id, edited)
        st.toast("Saved features")
