# src/abi_sauce/services/feature_ops.py
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

try:  # Biopython is optional at runtime in this project
    from Bio.SeqFeature import SeqFeature  # type: ignore[import-not-found]
except Exception:  # pragma: no cover

    class SeqFeature:  # type: ignore[no-redef]
        """Fallback placeholder when Biopython isn't installed."""


@dataclass
class Feature:
    """Typed feature with 1-based inclusive coordinates."""

    type: str
    start: int
    end: int
    strand: int = 0  # -1, 0, +1
    qualifiers: dict[str, Any] = field(default_factory=dict)
    location_text: str = ""


def _qnorm(q: dict[str, Any] | None) -> dict[str, list[Any]]:
    """Normalize qualifier dict so every value is a list."""
    out: dict[str, list[Any]] = {}
    for k, v in (q or {}).items():
        if isinstance(v, list):
            out[k] = v
        else:
            out[k] = [v]
    return out


def feature_from_biopython(f: SeqFeature, seq_len: int) -> Feature:
    """
    Best-effort conversion from Bio.SeqFeature to our typed Feature.

    Biopython uses 0-based, end-exclusive coordinates.
    We convert to 1-based, end-inclusive.
    """
    loc_txt = ""
    strand = 0
    start1 = 1
    end1 = seq_len

    try:
        loc = getattr(f, "location", None)
        if loc is not None:
            loc_txt = str(loc)
            raw_strand = getattr(loc, "strand", 0) or 0
            strand = 1 if raw_strand > 0 else (-1 if raw_strand < 0 else 0)

            # Handle both simple and compound locations.
            starts: list[int] = []
            ends: list[int] = []
            if hasattr(loc, "parts") and loc.parts:
                for p in loc.parts:  # type: ignore[attr-defined]
                    try:
                        starts.append(int(getattr(p, "start", 0)))
                        ends.append(int(getattr(p, "end", seq_len)))
                    except Exception:  # pragma: no cover
                        continue
            else:
                starts.append(int(getattr(loc, "start", 0)))
                ends.append(int(getattr(loc, "end", seq_len)))

            if starts:
                start1 = min(starts) + 1  # 0-based -> 1-based
            if ends:
                end1 = max(ends)  # end-exclusive -> inclusive by leaving as-is
    except Exception:
        # Keep safe defaults on any parsing issue.
        pass

    return Feature(
        type=str(getattr(f, "type", "feature")),
        start=start1,
        end=end1,
        strand=strand,
        qualifiers=_qnorm(getattr(f, "qualifiers", {}) or {}),
        location_text=loc_txt,
    )


def feature_to_dict(f: Feature) -> dict[str, Any]:
    """Dict for UI/JSON (keeps location text for fidelity)."""
    return {
        "type": f.type,
        "location": f.location_text or f"{f.start}..{f.end}",
        "start": f.start,
        "end": f.end,
        "strand": f.strand,
        "qualifiers": f.qualifiers,
    }


def features_to_dicts(items: Iterable[Feature]) -> list[dict[str, Any]]:
    return [feature_to_dict(x) for x in items]


def dict_to_feature(d: dict[str, Any]) -> Feature:
    return Feature(
        type=str(d.get("type", "feature")),
        start=int(d.get("start", 1)),
        end=int(d.get("end", d.get("start", 1))),
        strand=int(d.get("strand", 0)),
        qualifiers=_qnorm(d.get("qualifiers") or {}),
        location_text=str(
            d.get("location") or f"{d.get('start', 1)}..{d.get('end', 1)}"
        ),
    )


def normalize_feature_dicts(items: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Normalize existing dict features (strand -> int, qualifiers -> list[Any]).
    """
    out: list[dict[str, Any]] = []
    for d in items or []:
        f = dict_to_feature(d)
        out.append(feature_to_dict(f))
    return out


def auto_label(f: Feature) -> str:
    """Pick a nice display label from common qualifier keys."""
    q = f.qualifiers
    for k in ("label", "gene", "product", "note"):
        v = q.get(k)
        if isinstance(v, list) and v:
            return str(v[0])
        if isinstance(v, str) and v:
            return v
    return f.type
