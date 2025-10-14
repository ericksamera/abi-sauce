#!/usr/bin/env python3
from __future__ import annotations
import streamlit as st

from abi_sauce.models import Sample
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager


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
