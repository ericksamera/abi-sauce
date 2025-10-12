from __future__ import annotations
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.ui.controls import asset_selector
from abi_sauce.ui.components import asset_detail


def viewer_page():
    st.title("👀 Viewer")
    st.caption(
        "Pick any uploaded file to preview. AB1 shows chromatogram; FASTA/GenBank/ApE show sequence/features."
    )

    fm: FileManager = st.session_state._manager

    selected_id = asset_selector(
        fm, label="File", kinds=("Sequences", "Traces"), sidebar=True
    )
    if not selected_id:
        return

    asset = fm.get(selected_id)
    asset_detail(asset)
