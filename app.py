from __future__ import annotations

import streamlit as st

from abi_sauce.services.batch import (
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

_UPLOADER_NONCE_SESSION_KEY = "abi_sauce.uploader_nonce"
_UPLOADER_KEY_PREFIX = "abi_sauce.active_batch_uploads"


st.set_page_config(page_title="ABI Sauce", layout="wide")


def _active_uploader_key() -> str:
    uploader_nonce = int(st.session_state.get(_UPLOADER_NONCE_SESSION_KEY, 0))
    return f"{_UPLOADER_KEY_PREFIX}.{uploader_nonce}"


def _clear_active_batch_and_uploader() -> None:
    clear_active_batch(st.session_state)
    st.session_state[_UPLOADER_NONCE_SESSION_KEY] = (
        int(st.session_state.get(_UPLOADER_NONCE_SESSION_KEY, 0)) + 1
    )


def _sync_active_batch_from_uploader() -> None:
    uploaded_files = st.session_state.get(_active_uploader_key())
    if not uploaded_files:
        clear_active_batch(st.session_state)
        return

    uploads = normalize_uploaded_files(uploaded_files)
    next_signature = build_batch_signature(uploads)
    active_parsed_batch = get_active_parsed_batch(st.session_state)
    if (
        active_parsed_batch is not None
        and next_signature == get_active_batch_signature(st.session_state)
    ):
        return

    with st.spinner("Parsing uploaded ABI files..."):
        set_active_parsed_batch(st.session_state, parse_uploads(uploads))


st.sidebar.header("Active batch")

if st.sidebar.button("Clear active batch"):
    _clear_active_batch_and_uploader()
    st.rerun()

st.sidebar.file_uploader(
    "Upload ABI trace files",
    type=["ab1", "abi"],
    accept_multiple_files=True,
    key=_active_uploader_key(),
    help="Loaded files persist across pages for this session.",
)

_sync_active_batch_from_uploader()

active_parsed_batch = get_active_parsed_batch(st.session_state)
if active_parsed_batch is None:
    st.sidebar.info("No active batch loaded.")
else:
    st.sidebar.write(
        {
            "uploaded_files": len(active_parsed_batch.uploads),
            "parsed_files": len(active_parsed_batch.parsed_records),
            "failed_files": len(active_parsed_batch.parse_errors),
        }
    )
    if active_parsed_batch.parse_errors:
        st.sidebar.warning(
            f"{len(active_parsed_batch.parse_errors)} file(s) failed to parse."
        )

navigation = st.navigation(
    [
        st.Page("pages/00_upload_and_parse.py", title="Upload and Parse"),
        st.Page("pages/01_chromatogram_preview.py", title="Chromatogram Viewer"),
    ]
)

navigation.run()
