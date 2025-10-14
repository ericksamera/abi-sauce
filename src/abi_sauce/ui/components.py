# src/abi_sauce/ui/components.py
from __future__ import annotations
from typing import List, Optional
import streamlit as st

from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset
from abi_sauce.services.sample_manager import SampleManager

# New focused modules
from abi_sauce.ui.trace_viewer import render_trace_asset
from abi_sauce.ui.sequence_viewer import render_sequence_asset
from abi_sauce.ui.sample_editor import sample_editor  # re-exported


def asset_table(assets: List[AssetBase]) -> Optional[str]:
    if not assets:
        st.info("No files uploaded yet.")
        return None
    labels = [f"[{a.kind.value}] {a.name} • {getattr(a, 'length', '')}" for a in assets]
    ids = [a.id for a in assets]
    selected = st.radio("Assets", ids, format_func=lambda x: labels[ids.index(x)])
    return selected


def asset_detail(asset: AssetBase) -> None:
    """
    Delegates to the appropriate, focused UI renderer.
    This function name remains unchanged for back-compat with viewer_page.
    """
    if isinstance(asset, SequenceAsset):
        render_sequence_asset(asset)
        return
    if isinstance(asset, TraceAsset):
        render_trace_asset(asset)
        return
    # Fallback for unknown asset types
    st.write(asset)


def sample_table(samples: List["SampleManager"]) -> Optional[str]:
    # Signature retained; behavior unchanged.
    if not samples:
        st.info("No samples yet — upload files on the Uploads page.")
    labels = [f"{s.name}" for s in samples]
    ids = [s.id for s in samples]
    selected = (
        st.radio("Samples", ids, format_func=lambda x: labels[ids.index(x)])
        if samples
        else None
    )
    return selected


__all__ = [
    "asset_table",
    "asset_detail",
    "sample_table",
    "sample_editor",  # from ui.sample_editor
]
