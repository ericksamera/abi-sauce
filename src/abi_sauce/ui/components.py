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
    Renders the detail view for an asset.
    - SequenceAsset: shows sequence (+ features) and FASTA.
    - TraceAsset: chromatogram plot (via build_trace_fig), Mott trim preview,
      download trimmed FASTA, and **Apply trim to Sample**.
    Looks up SampleManager from st.session_state so Viewer/Uploads work unchanged.
    """
    # --- Sequence assets ---
    if isinstance(asset, SequenceAsset):
        st.subheader("Sequence")
        st.code(
            asset.sequence[:5000] + ("…" if getattr(asset, "length", 0) > 5000 else "")
        )
        with st.expander("Features"):
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

    # --- Trace assets (AB1) ---
    st.subheader("Trace (AB1)")
    seq = asset.sequence or ""
    ploc = asset.base_positions or []
    quals = asset.qualities or []
    channels = {k: (asset.channels.get(k) or []) for k in ["A", "C", "G", "T"]}

    # Called sequence (if present)
    if seq:
        st.text(f"Sequence length: {len(seq)} bases")
        fa = asset.to_fasta()
        if fa:
            st.download_button(
                "Download called sequence (FASTA)",
                data=fa,
                file_name=f"{asset.name}.fasta",
            )

    if not any(channels.values()):
        st.info("No channel data found.")
        return

    # --- Plot options ---
    with st.expander("Chart options", expanded=False):
        mode = st.radio("X-axis mode", ["samples", "bases"], index=0, horizontal=True)
        show_rangeslider = st.toggle("Show range slider", value=True)
        show_grid = st.toggle("Show grid", value=True)
        peak_only_hover = st.toggle("Hover only on peaks", value=True)

    # --- Trim settings (Mott) ---
    with st.expander("Trim settings (Mott)", expanded=True):
        err_key = f"mott_err_{asset.id}"
        minlen_key = f"mott_minlen_{asset.id}"
        err_limit = st.number_input(
            "Error limit (typ. 0.05)", 0.0, 0.5, 0.05, 0.01, key=err_key
        )
        min_len = st.number_input("Minimum kept length", 1, 1000, 20, 1, key=minlen_key)

    # Compute trim (0-based inclusive l,r) if possible
    l_idx = r_idx = None
    trimmed_seq = None
    if seq and quals:
        out = trim_trace_asset_mott(seq, quals, error_limit=err_limit, min_len=min_len)
        if out:
            l_idx, r_idx, trimmed_seq = out

    # Build figure
    seq_len = len(seq) if seq else len(ploc)
    trim_range = None  # 1-based inclusive for build_trace_fig
    if l_idx is not None and r_idx is not None and l_idx <= r_idx:
        trim_range = (int(l_idx) + 1, int(r_idx) + 1)

    fig = build_trace_fig(
        channels=channels,
        ploc=ploc,
        seq_len=seq_len,
        quals=quals,
        seq=seq or None,
        mode=mode,
        resample_method="max",
        show_trim=True,
        trim_range=trim_range,
        uirevision=f"{asset.id}_chrom_v2",
        height=360,
        show_grid=show_grid,
        peak_only_hover=peak_only_hover,
        show_rangeslider=show_rangeslider,
    )
    st.plotly_chart(fig, use_container_width=True)

    # --- Download trimmed FASTA (if available) ---
    if trimmed_seq:
        fasta = f">{asset.name} | Mott {l_idx+1}-{r_idx+1}\n"
        fasta += "\n".join(
            trimmed_seq[i : i + 70] for i in range(0, len(trimmed_seq), 70)
        )
        st.download_button(
            "Download Mott-trimmed FASTA",
            data=fasta,
            file_name=f"{asset.name}.mott.fasta",
        )

    # --- Apply Mott trim to Sample (if there are samples that reference this asset) ---
    sm: Optional[SampleManager] = st.session_state.get(
        "_samples", None
    )  # auto-discovered
    if sm and trimmed_seq:
        samples_for_asset = [s for s in sm.list() if asset.id in (s.asset_ids or [])]

        if not samples_for_asset:
            st.info(
                "No Sample references this file yet. Add it to a Sample on the Samples page."
            )
            return

        with st.container():
            st.markdown("#### Apply to Sample")
            if len(samples_for_asset) == 1:
                target = samples_for_asset[0]

                def _apply_one():
                    sm.set_sequence_override(target.id, trimmed_seq)
                    st.toast(f"Applied trim to “{target.name}”")

                st.button(f"Apply Mott trim to “{target.name}”", on_click=_apply_one)
            else:
                sel_key = f"trim_target_{asset.id}"
                # Default to the first one
                default_idx = 0
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
                    index=default_idx,
                    format_func=lambda sid: names[ids.index(sid)],
                    key=sel_key,
                )
                st.button("Apply Mott trim", on_click=_apply_many)


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
        # default radio value
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
