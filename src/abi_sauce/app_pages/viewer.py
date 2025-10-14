from __future__ import annotations

import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.ui.components import asset_detail
from abi_sauce.ui.controls import asset_selector


def viewer_page():

    fm: FileManager = st.session_state._manager

    selected_id = asset_selector(
        fm, label="File", kinds=("Sequences", "Traces"), sidebar=True
    )
    if not selected_id:
        return

    asset = fm.get(selected_id)
    asset_detail(asset)
