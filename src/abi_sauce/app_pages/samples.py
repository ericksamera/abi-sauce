from __future__ import annotations
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.ui.controls import sample_selector
from abi_sauce.ui.components import sample_editor


def samples_page():
    st.title("Samples")
    st.caption(
        "Each upload becomes an editable Sample. Edits never mutate the original file."
    )

    fm: FileManager = st.session_state._manager
    sm: SampleManager = st.session_state._samples

    sid = sample_selector(sm, label="Samples")
    if not sid:
        return

    sample = sm.get(sid)
    sample_editor(sample, fm, sm)
