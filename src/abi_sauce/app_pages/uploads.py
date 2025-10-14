# src/abi_sauce/app_pages/uploads.py
from __future__ import annotations

import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager


def uploads_page() -> None:
    st.title("Uploads")

    fm: FileManager = st.session_state._manager
    sm: SampleManager = st.session_state._samples

    st.caption(
        "Drop FASTA / GenBank / ApE / AB1 files. "
        "Each new asset automatically becomes a Sample for editing later."
    )

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
            except Exception as e:  # narrow when you add specific errors
                st.warning(f"{f.name}: {e}")

    with st.sidebar:
        st.header("Actions")

        def _clear_all() -> None:
            fm.clear()
            st.session_state["_assets_cleared"] = True
            st.rerun()

        st.button("Clear all", on_click=_clear_all, use_container_width=True)
