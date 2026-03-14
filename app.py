from __future__ import annotations

from collections.abc import Sequence
from urllib.parse import quote

import streamlit as st

from abi_sauce.batch_download_ui import render_batch_download_controls
from abi_sauce.viewer_state import get_batch_trim_state
from abi_sauce.services.batch import apply_trim_configs, prepare_batch_download
from abi_sauce.trim_state import resolve_batch_trim_inputs

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
    get_active_uploads,
    merge_uploads,
    set_active_parsed_batch,
)
from abi_sauce.viewer_state import clear_viewer_session_state
from abi_sauce.assembly_state import clear_assembly_session_state

_UPLOADER_NONCE_SESSION_KEY = "abi_sauce.uploader_nonce"
_UPLOADER_KEY_PREFIX = "abi_sauce.active_batch_uploads"
_SIDEBAR_DOWNLOAD_POPOVER_KEY = "abi_sauce.sidebar.download_popover"


st.set_page_config(page_title="ABI Sauce", layout="wide")
st.markdown(
    """
<style>
section.stMain .block-container {
    padding-top: 2.5rem;
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
    clear_assembly_session_state(st.session_state)
    _bump_uploader_nonce()


@st.dialog(
    "Restart session",
    width="small",
    on_dismiss="rerun",
)
def _restart_session_dialog() -> None:
    st.error(
        "This will clear the active batch, selected sample, and shared trim/viewer state."
    )

    confirm_col, cancel_col = st.columns(2)

    with confirm_col:
        confirm_clicked = st.button(
            "Confirm restart",
            type="primary",
            use_container_width=True,
        )

    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
        )

    if cancel_clicked:
        st.rerun()

    if confirm_clicked:
        _clear_active_batch_and_uploader()
        st.rerun()


def _activate_uploaded_files(uploaded_files: Sequence[UploadedFileLike]) -> None:
    uploads = merge_uploads(
        get_active_uploads(st.session_state),
        normalize_uploaded_files(uploaded_files),
    )
    next_signature = build_batch_signature(uploads)
    if next_signature == get_active_batch_signature(st.session_state):
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
        "Load ABI trace files into the shared session batch. New filenames are appended; reloading the same filename replaces that entry."
    )
    uploaded_files = st.file_uploader(
        "Choose one or more ABI trace files",
        type=["ab1", "abi"],
        accept_multiple_files=True,
        key=_active_uploader_key(),
    )
    # if uploaded_files:
    #     st.write(
    #         {
    #             "selected_files": len(uploaded_files),
    #             "filenames": [uploaded_file.name for uploaded_file in uploaded_files],
    #         }
    #     )

    load_col, cancel_col = st.columns(2)
    with load_col:
        load_clicked = st.button(
            "Add to session",
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


def _prepare_active_batch():
    parsed_batch = get_active_parsed_batch(st.session_state)
    if parsed_batch is None or not parsed_batch.parsed_records:
        return None

    trim_state = get_batch_trim_state(st.session_state)
    resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)
    return apply_trim_configs(
        parsed_batch,
        default_trim_config=resolved_trim_inputs.default_trim_config,
        trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
    )


def _build_blast_url() -> str | None:
    prepared_batch = _prepare_active_batch()
    if prepared_batch is None:
        return None

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fasta",
        concatenate_batch=True,
        filename_stem="abi-sauce-batch",
        require_min_length=True,
        fasta_line_width=None,
    )
    if not artifact.is_downloadable or not isinstance(artifact.data, str):
        return None

    query_text = artifact.data.strip()
    return (
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        "?PROGRAM=blastn"
        "&PAGE_TYPE=BlastSearch"
        "&LINK_LOC=blasthome"
        f"&QUERY={quote(query_text, safe='')}"
    )


def _render_sidebar_download_popover() -> None:
    prepared_batch = _prepare_active_batch()
    if prepared_batch is None:
        return

    with st.popover(
        "Export Sequences",
        key=_SIDEBAR_DOWNLOAD_POPOVER_KEY,
        width="stretch",
        icon=":material/download:",
    ):
        render_batch_download_controls(
            prepared_batch=prepared_batch,
            key_prefix="abi_sauce.sidebar_download",
            default_filename_stem="abi-sauce-trim",
            button_label="Export Sequences",
            compact=True,
        )


active_parsed_batch = get_active_parsed_batch(st.session_state)

pages = [
    st.Page("pages/00_home.py", title="Home", icon=":material/home:"),
]

if active_parsed_batch is not None and active_parsed_batch.parsed_records:
    pages.extend(
        [
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
            st.Page(
                "pages/02_reference_alignment.py",
                title="Reference Alignment",
                icon=":material/compare_arrows:",
            ),
            st.Page(
                "pages/03_assembly.py",
                title="Assembly",
                icon=":material/merge_type:",
            ),
        ]
    )

current_page = st.navigation(
    pages,
    position="top",
)

st.sidebar.title(":apple: abi-sauce")
st.sidebar.caption(
    "Use the shared active batch across Home, Sample Viewer, Batch Viewer, and Assembly."
)

load_button_label = "Upload Files" if active_parsed_batch is None else "Upload Files"
with st.sidebar:
    col_import, col_reset = st.columns([4, 1], border=False)

    with col_import:
        if st.button(load_button_label, width="stretch"):
            _upload_batch_dialog()

    with col_reset:
        if st.button(
            "",
            type="primary",
            key="clear_batch",
            icon=":material/delete_sweep:",
            width="stretch",
            disabled=active_parsed_batch is None,
        ):
            _restart_session_dialog()

    is_assembly_page = current_page.title == "Assembly"
    if (
        active_parsed_batch is not None
        and active_parsed_batch.parsed_records
        and not is_assembly_page
    ):
        _render_sidebar_download_popover()

        st.divider()

        blast_url = _build_blast_url()
        if blast_url is None:
            st.caption(
                "No trimmed sequences meeting the minimum-length requirement are available for BLAST."
            )
        elif len(blast_url) > 60_000:
            st.caption(
                "Batch is too large for a reliable BLAST URL; export FASTA and paste/upload it in BLAST instead."
            )
        else:
            st.link_button(
                "BLAST Trimmed Sequences",
                blast_url,
                width="stretch",
                icon=":material/rocket_launch:",
            )

current_page.run()
