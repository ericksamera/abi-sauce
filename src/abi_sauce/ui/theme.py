# src/abi_sauce/ui/theme.py
from __future__ import annotations

from typing import Final

# ---- Colors ----
BASE_COLORS: Final[dict[str, str]] = {
    "A": "#16a34a",  # green-600
    "C": "#2563eb",  # blue-600
    "G": "#111827",  # gray-900
    "T": "#dc2626",  # red-600
}
DEFAULT_TEXT_COLOR: Final[str] = "#111827"
MUTED_TEXT_COLOR: Final[str] = "#666666"
UNKNOWN_COLOR: Final[str] = "#888888"

# ---- Hover / markers ----
HOVER_SPIKE_COLOR: Final[str] = "#444444"
HOVER_SPIKE_THICKNESS: Final[int] = 2
RANGE_SLIDER_THICKNESS: Final[float] = 0.12
WAVE_OPACITY: Final[float] = 0.95
QUALITY_BAR_OPACITY: Final[float] = 0.18

# ---- Hatched shading (trim flanks) ----
HATCH_BASE_ALPHA: Final[float] = 0.15
HATCH_STRIPE_ALPHA: Final[float] = 0.30

# ---- Layout defaults ----
MARGIN_DEFAULT: Final[dict[str, int]] = {"l": 40, "r": 40, "t": 40, "b": 40}
LEGEND_DEFAULT: Final[dict[str, object]] = {
    "orientation": "h",
    "yanchor": "bottom",
    "y": 1.02,
    "xanchor": "left",
    "x": 0,
}


def base_color(base: str) -> str:
    """Color for a given base letter (fallback to default text color)."""
    return BASE_COLORS.get(base.upper(), DEFAULT_TEXT_COLOR)


def letter_font_size_by_columns(ncols: int) -> int:
    """Keep per-column letters readable on dense alignments."""
    if ncols <= 600:
        return 9
    if ncols <= 1200:
        return 8
    if ncols <= 1800:
        return 7
    return 6
