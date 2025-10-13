# src/abi_sauce/ui/defaults.py
from __future__ import annotations
from typing import Dict, Any
import streamlit as st

# ---- Shared UI defaults ----
CHART_HEIGHT_DEFAULT = 560
INITIAL_BASE_SPAN = 40
BASE_TICK_EVERY = 10
TRIM_ERROR_DEFAULT = 0.05
TRIM_MINLEN_DEFAULT = 20

# Keys are suffixes appended to a per-asset prefix, e.g. f"{prefix}_rangeslider"
DEFAULT_WIDGETS: Dict[str, Any] = {
    "rangeslider": True,
    "grid": True,
    "hover_peaks_only": True,
    "err": TRIM_ERROR_DEFAULT,
    "minlen": TRIM_MINLEN_DEFAULT,
    "height": CHART_HEIGHT_DEFAULT,
    "initialized": False,  # used to apply initial x-range once
}


def init_state(prefix: str, overrides: Dict[str, Any] | None = None) -> None:
    """
    Ensure all default widget values exist in st.session_state
    for a given prefix (e.g., per-asset).
    """
    vals = {**DEFAULT_WIDGETS, **(overrides or {})}
    for k, v in vals.items():
        st.session_state.setdefault(f"{prefix}_{k}", v)
