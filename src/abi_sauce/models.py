# src/abi_sauce/models.py
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any
from datetime import datetime, timezone


class AssetKind(str, Enum):
    SEQUENCE = "sequence"  # e.g., FASTA/GenBank/ApE
    TRACE = "trace"  # e.g., AB1 (Sanger chromatogram)


@dataclass
class AssetBase:
    id: str
    name: str
    ext: str
    size: int
    checksum_md5: str
    created_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc), init=False
    )
    kind: AssetKind = field(init=False)


@dataclass
class SequenceAsset(AssetBase):
    sequence: str  # plain sequence (A/C/G/T, possibly N)
    description: str = ""
    length: int = 0
    meta: Dict[str, Any] = field(default_factory=dict)  # e.g., source, locus, etc.
    features: List[Dict[str, Any]] = field(
        default_factory=list
    )  # simplified feature dicts

    def __post_init__(self):
        self.kind = AssetKind.SEQUENCE
        if not self.length:
            self.length = len(self.sequence)

    def to_fasta(self, header: Optional[str] = None) -> str:
        hdr = header or (self.description or self.name)
        seq = self.sequence
        # wrap 70-chars per FASTA convention
        wrapped = "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))
        return f">{hdr}\n{wrapped}\n"


@dataclass
class TraceAsset(AssetBase):
    sequence: Optional[str] = None
    qualities: Optional[List[int]] = None  # per-base PHRED, if present
    base_positions: Optional[List[int]] = None  # peak indices (PLOC2), if present
    channels: Dict[str, Optional[List[int]]] = field(
        default_factory=lambda: {"A": None, "C": None, "G": None, "T": None}
    )
    meta: Dict[str, Any] = field(
        default_factory=dict
    )  # e.g., instrument, run date, ABI tags subset

    def __post_init__(self):
        self.kind = AssetKind.TRACE

    def to_fasta(self, header: Optional[str] = None) -> Optional[str]:
        if not self.sequence:
            return None
        hdr = header or self.name
        wrapped = "\n".join(
            self.sequence[i : i + 70] for i in range(0, len(self.sequence), 70)
        )
        return f">{hdr}\n{wrapped}\n"


@dataclass
class Sample:
    id: str
    name: str
    asset_ids: List[str] = field(default_factory=list)
    description: str = ""
    tags: List[str] = field(default_factory=list)
    primary_asset_id: Optional[str] = None

    # user edits (do not mutate original assets)
    sequence_override: Optional[str] = None
    feature_overrides: Optional[List[Dict[str, Any]]] = None

    created_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc), init=False
    )

    # helpers (these do not persist assets; managers should pass in assets as a dict)
    def effective_sequence(self, assets: Dict[str, AssetBase]) -> Optional[str]:
        if self.sequence_override:
            return self.sequence_override
        # choose primary, else first sequence-like thing with bases
        cand_id = self.primary_asset_id or next(
            (aid for aid in self.asset_ids if aid in assets), None
        )
        if not cand_id:
            return None
        a = assets.get(cand_id)
        if isinstance(a, SequenceAsset):
            return a.sequence
        if isinstance(a, TraceAsset):
            return a.sequence
        return None

    def effective_features(self, assets: Dict[str, AssetBase]) -> List[Dict[str, Any]]:
        if self.feature_overrides is not None:
            return self.feature_overrides
        cand_id = self.primary_asset_id or next(
            (aid for aid in self.asset_ids if aid in assets), None
        )
        a = assets.get(cand_id)
        if isinstance(a, SequenceAsset) and a.features:
            return a.features
        return []

    def to_fasta(self, assets: Dict[str, AssetBase]) -> Optional[str]:
        seq = self.effective_sequence(assets)
        if not seq:
            return None
        wrapped = "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))
        return f">{self.name}\n{wrapped}\n"
