# src/abi_sauce/app_pages/uploads.py
from __future__ import annotations
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.ui.components import asset_table, asset_detail

def uploads_page():
    st.title("📥 Uploads")

    fm: FileManager = st.session_state._manager
    sm: SampleManager = st.session_state._samples

    st.caption("Drop FASTA / GenBank / ApE / AB1 files. Each new asset automatically becomes a Sample for editing later.")

    uploaded = st.file_uploader(
        "Drop files here",
        type=["fa", "fasta", "fna", "gb", "gbk", "gbff", "ape", "ab1", "abi"],
        accept_multiple_files=True,
    )
    if uploaded:
        for f in uploaded:
            raw = f.getvalue()
            try:
                sample_ids = sm.import_bytes(name=f.name, raw=raw)
                if sample_ids:
                    st.toast(f"Imported {f.name} → {len(sample_ids)} sample(s)")
                else:
                    st.toast(f"Skipped duplicate {f.name}")
            except Exception as e:
                st.warning(f"{f.name}: {e}")

    st.subheader("Assets")
    assets = fm.list()
    selected_id = asset_table(assets)
    if selected_id:
        asset = fm.get(selected_id)
        asset_detail(asset)

    with st.sidebar:
        st.header("Actions")

        def _clear_all():
            fm.clear()
            # set a small session_state flag so UI reflects cleared state on rerun
            st.session_state["_assets_cleared"] = True

        st.button("Clear all", width="stretch", on_click=_clear_all)
