from __future__ import annotations

from collections.abc import Sequence

import streamlit as st

from abi_sauce.services.batch import (
    UploadedFileLike,
    build_batch_signature,
    normalize_uploaded_files,
    parse_uploads,
)
from abi_sauce.upload_state import (
    clear_active_batch,
    get_active_batch_signature,
    get_active_parsed_batch,
    set_active_parsed_batch,
)
from abi_sauce.viewer_state import clear_viewer_session_state

_UPLOADER_NONCE_SESSION_KEY = "abi_sauce.uploader_nonce"
_UPLOADER_KEY_PREFIX = "abi_sauce.active_batch_uploads"


st.set_page_config(page_title="ABI Sauce", layout="wide")
st.markdown(
    """
<style>
section.stMain .block-container {
    padding-top: 2rem;
    z-index: 1;
}
</style>""",
    unsafe_allow_html=True,
)


def _active_uploader_key() -> str:
    uploader_nonce = int(st.session_state.get(_UPLOADER_NONCE_SESSION_KEY, 0))
    return f"{_UPLOADER_KEY_PREFIX}.{uploader_nonce}"


def _bump_uploader_nonce() -> None:
    st.session_state[_UPLOADER_NONCE_SESSION_KEY] = (
        int(st.session_state.get(_UPLOADER_NONCE_SESSION_KEY, 0)) + 1
    )


def _clear_active_batch_and_uploader() -> None:
    clear_active_batch(st.session_state)
    clear_viewer_session_state(st.session_state)
    _bump_uploader_nonce()


def _activate_uploaded_files(uploaded_files: Sequence[UploadedFileLike]) -> None:
    uploads = normalize_uploaded_files(uploaded_files)
    next_signature = build_batch_signature(uploads)
    active_parsed_batch = get_active_parsed_batch(st.session_state)
    if (
        active_parsed_batch is not None
        and next_signature == get_active_batch_signature(st.session_state)
    ):
        _bump_uploader_nonce()
        return

    with st.spinner("Parsing uploaded ABI files..."):
        set_active_parsed_batch(st.session_state, parse_uploads(uploads))
    _bump_uploader_nonce()


@st.dialog(
    "Load ABI trace files",
    width="small",
    on_dismiss="rerun",
)
def _upload_batch_dialog() -> None:
    st.caption(
        "Load a shared active batch for the Home, Sample Viewer, and Batch Viewer pages."
    )
    uploaded_files = st.file_uploader(
        "Choose one or more ABI trace files",
        type=["ab1", "abi"],
        accept_multiple_files=True,
        key=_active_uploader_key(),
    )
    if uploaded_files:
        st.write(
            {
                "selected_files": len(uploaded_files),
                "filenames": [uploaded_file.name for uploaded_file in uploaded_files],
            }
        )

    load_col, cancel_col = st.columns(2)
    with load_col:
        load_clicked = st.button(
            "Load batch",
            type="primary",
            use_container_width=True,
            disabled=not uploaded_files,
        )
    with cancel_col:
        cancel_clicked = st.button("Cancel", use_container_width=True)

    if cancel_clicked:
        _bump_uploader_nonce()
        st.rerun()

    if load_clicked:
        if not uploaded_files:
            st.warning("Choose one or more .ab1 files to load.")
            return
        _activate_uploaded_files(uploaded_files)
        st.rerun()


active_parsed_batch = get_active_parsed_batch(st.session_state)

st.sidebar.title("Workspace")
st.sidebar.caption(
    "Use the shared active batch across Home, Sample Viewer, and Batch Viewer."
)

load_button_label = (
    "Load ABI files" if active_parsed_batch is None else "Replace active batch"
)
if st.sidebar.button(load_button_label, type="primary", use_container_width=True):
    _upload_batch_dialog()

if st.sidebar.button(
    "Clear active batch",
    use_container_width=True,
    disabled=active_parsed_batch is None,
):
    _clear_active_batch_and_uploader()
    st.rerun()

navigation = st.navigation(
    [
        st.Page("pages/00_home.py", title="Home", icon=":material/home:"),
        st.Page(
            "pages/01_chromatogram_preview.py",
            title="Sample Viewer",
            icon=":material/biotech:",
        ),
        st.Page(
            "pages/00_upload_and_parse.py",
            title="Batch Viewer",
            icon=":material/view_list:",
        ),
    ],
    position="top",
)

navigation.run()
