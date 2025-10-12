# src/abi_sauce/ui/controls.py
from __future__ import annotations
from typing import Iterable, Optional, Sequence
import streamlit as st

from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager


def asset_selector(
    fm: FileManager,
    *,
    label: str = "File",
    kinds: Sequence[str] = ("Sequences", "Traces"),
    default_id: Optional[str] = None,
    sidebar: bool = True,
) -> Optional[str]:
    """Unified asset picker with type filter + search."""
    assets = fm.list()
    if not assets:
        st.info("No files uploaded yet.")
        return None

    target = st.sidebar if sidebar else st

    with target:
        target.subheader("Select file")
        type_choices = ["Sequences", "Traces"]
        enabled = target.multiselect("Types", type_choices, default=list(kinds))
        query = target.text_input("Search by name", "")

    def _match(a: AssetBase) -> bool:
        kind_ok = (isinstance(a, SequenceAsset) and "Sequences" in enabled) or (
            isinstance(a, TraceAsset) and "Traces" in enabled
        )
        q_ok = (query.lower() in a.name.lower()) if query else True
        return kind_ok and q_ok

    filtered = [a for a in assets if _match(a)]
    if not filtered:
        st.info("No files match your filters.")
        return None

    ids = [a.id for a in filtered]

    def _label(a: AssetBase) -> str:
        if isinstance(a, SequenceAsset):
            extra = f"{a.length} bp"
            kind = "sequence"
        elif isinstance(a, TraceAsset):
            extra = f"{len(a.sequence) if a.sequence else 0} bases"
            kind = "trace"
        else:
            extra, kind = "", "asset"
        return f"[{kind}] {a.name} • {extra}".strip()

    labels = {a.id: _label(a) for a in filtered}

    default_index = 0
    if default_id and default_id in ids:
        default_index = ids.index(default_id)
    elif (
        "_viewer_selected" in st.session_state
        and st.session_state._viewer_selected in ids
    ):
        default_index = ids.index(st.session_state._viewer_selected)

    picked = target.selectbox(
        label, ids, index=default_index, format_func=lambda _id: labels[_id]
    )
    st.session_state._viewer_selected = picked
    return picked


def sample_selector(
    sm: SampleManager,
    *,
    label: str = "Sample",
    exclude: Iterable[str] = (),
    default_id: Optional[str] = None,
    sidebar: bool = False,
) -> Optional[str]:
    """Unified sample picker."""
    samples = [s for s in sm.list() if s.id not in set(exclude)]
    if not samples:
        st.info("No samples yet — upload files on the Uploads page.")
        return None

    target = st.sidebar if sidebar else st
    ids = [s.id for s in samples]
    labels = {s.id: s.name for s in samples}

    idx = 0
    if default_id and default_id in ids:
        idx = ids.index(default_id)

    return target.radio(label, ids, index=idx, format_func=lambda sid: labels[sid])
