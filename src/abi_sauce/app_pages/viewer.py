from __future__ import annotations
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.ui.components import asset_detail
from abi_sauce.models import SequenceAsset, TraceAsset, AssetBase


def viewer_page():
    st.title("👀 Viewer")
    st.caption(
        "Pick any uploaded file from the sidebar to preview it. AB1 shows chromatogram + sequence; FASTA/GenBank/ApE show sequence (and features if present)."
    )

    fm: FileManager = st.session_state._manager
    assets = fm.list()
    if not assets:
        st.info(
            "No files yet — go to the Uploads page and add FASTA / GenBank / ApE / AB1."
        )
        return

    # Sidebar selector + filters
    with st.sidebar:
        st.header("Select file")

        # Type filter
        type_choices = ["Sequences", "Traces"]
        selected_types = st.multiselect(
            "Types",
            type_choices,
            default=type_choices,  # show everything by default
        )

        # Quick search
        q = st.text_input("Search by name", "")

        def _include(a: AssetBase) -> bool:
            is_seq = isinstance(a, SequenceAsset)
            is_trace = isinstance(a, TraceAsset)
            type_ok = (is_seq and "Sequences" in selected_types) or (
                is_trace and "Traces" in selected_types
            )
            q_ok = (q.lower() in a.name.lower()) if q else True
            return type_ok and q_ok

        filtered = [a for a in assets if _include(a)]
        if not filtered:
            st.info("No files match your filters.")
            return

        # Labels
        def _label(a: AssetBase) -> str:
            if isinstance(a, SequenceAsset):
                extra = f"{a.length} bp"
                kind = "sequence"
            elif isinstance(a, TraceAsset):
                extra = f"{len(a.sequence) if a.sequence else 0} bases"
                kind = "trace"
            else:
                extra = ""
                kind = "asset"
            return f"[{kind}] {a.name} • {extra}".strip()

        ids = [a.id for a in filtered]
        labels = {a.id: _label(a) for a in filtered}

        default_index = 0
        if (
            "_viewer_selected" in st.session_state
            and st.session_state._viewer_selected in ids
        ):
            default_index = ids.index(st.session_state._viewer_selected)

        selected_id = st.selectbox(
            "File",
            ids,
            index=default_index,
            format_func=lambda _id: labels[_id],
        )
        st.session_state._viewer_selected = selected_id

    # Main view
    asset = fm.get(selected_id)
    asset_detail(asset)
