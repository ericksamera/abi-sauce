from __future__ import annotations
import streamlit as st
from datetime import datetime

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.services.workspace import WorkspaceManager


def projects_page():
    st.title("Projects / Workspaces")

    fm: FileManager = st.session_state._manager
    sm: SampleManager = st.session_state._samples
    wm = WorkspaceManager()

    st.subheader("Export")
    snap = wm.export_json(fm, sm, indent=2)
    fname = f"abi-sauce-{datetime.now().strftime('%Y%m%d-%H%M%S')}.json"
    st.download_button(
        "Download current workspace",
        data=snap,
        file_name=fname,
        use_container_width=True,
    )

    st.divider()
    st.subheader("Import")
    clear_first = st.toggle("Clear existing items before import", value=True)
    upl = st.file_uploader("Choose a workspace JSON", type=["json"])
    if upl:
        try:
            payload = upl.getvalue().decode("utf-8")
            n_assets, n_samples = wm.import_json(
                fm, sm, payload, clear_first=clear_first
            )
            st.success(f"Imported {n_assets} assets and {n_samples} samples.")
            st.toast("Workspace imported")
        except Exception as e:
            st.error(f"Import failed: {e}")
