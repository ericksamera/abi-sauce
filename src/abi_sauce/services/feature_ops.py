# src/abi_sauce/services/feature_ops.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List

try:
    from Bio.SeqFeature import SeqFeature
except Exception:
    SeqFeature = object  # type: ignore


@dataclass
class Feature:
    """Typed feature (1-based inclusive coords)."""

    type: str
    start: int
    end: int
    strand: int = 0  # -1, 0, +1
    qualifiers: Dict[str, Any] = field(default_factory=dict)
    location_text: str = ""


def _qnorm(q: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in (q or {}).items():
        if isinstance(v, list):
            out[k] = v
        else:
            out[k] = [v]
    return out


def feature_from_biopython(f: "SeqFeature", seq_len: int) -> Feature:
    """Best-effort conversion from Bio.SeqFeature to our typed Feature."""
    loc_txt = ""
    strand = 0
    start1, end1 = 1, seq_len
    try:
        loc_txt = str(f.location)
        strand = int(getattr(f.location, "strand", 0) or 0)
        # Bio uses 0-based start, end-exclusive; convert to 1-based inclusive
        start1 = int(getattr(f.location, "start", 0)) + 1
        end1 = int(getattr(f.location, "end", seq_len))
    except Exception:
        pass
    return Feature(
        type=getattr(f, "type", "feature"),
        start=start1,
        end=end1,
        strand=strand,
        qualifiers=_qnorm(getattr(f, "qualifiers", {}) or {}),
        location_text=loc_txt,
    )


def feature_to_dict(f: Feature) -> Dict[str, Any]:
    """Dict for UI / JSON (keeps location text for fidelity)."""
    return {
        "type": f.type,
        "location": f.location_text or f"{f.start}..{f.end}",
        "start": f.start,
        "end": f.end,
        "strand": f.strand,
        "qualifiers": f.qualifiers,
    }


def features_to_dicts(items: Iterable[Feature]) -> List[Dict[str, Any]]:
    return [feature_to_dict(x) for x in items]


def dict_to_feature(d: Dict[str, Any]) -> Feature:
    return Feature(
        type=str(d.get("type", "feature")),
        start=int(d.get("start", 1)),
        end=int(d.get("end", d.get("start", 1))),
        strand=int(d.get("strand", 0)),
        qualifiers=_qnorm(d.get("qualifiers") or {}),
        location_text=str(d.get("location") or f"{d.get('start',1)}..{d.get('end',1)}"),
    )


def normalize_feature_dicts(items: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Normalize existing dict features (strand -> int, qualifiers -> list[str])."""
    out: List[Dict[str, Any]] = []
    for d in items or []:
        f = dict_to_feature(d)
        out.append(feature_to_dict(f))
    return out


def auto_label(f: Feature) -> str:
    """Pick a nice display label."""
    q = f.qualifiers
    for k in ("label", "gene", "product", "note"):
        v = q.get(k)
        if isinstance(v, list) and v:
            return str(v[0])
        if isinstance(v, str) and v:
            return v
    return f.type
