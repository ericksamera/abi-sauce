# src/abi_sauce/models.py
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import Enum
from typing import Any

from abi_sauce.types import BaseChar


class AssetKind(str, Enum):
    SEQUENCE = "sequence"
    TRACE = "trace"


@dataclass
class AssetBase:
    id: str
    name: str
    ext: str
    size: int
    checksum_md5: str
    created_at: datetime = field(default_factory=lambda: datetime.now(UTC))
    updated_at: datetime = field(default_factory=lambda: datetime.now(UTC))
    kind: AssetKind = field(init=False)

    def touch(self) -> None:
        self.updated_at = datetime.now(UTC)


@dataclass
class SequenceAsset(AssetBase):
    sequence: str
    description: str = ""
    length: int = 0
    meta: dict[str, Any] = field(default_factory=dict)
    features: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.kind = AssetKind.SEQUENCE
        if not self.length:
            self.length = len(self.sequence)

    def to_fasta(self, header: str | None = None) -> str:
        hdr = header or (self.description or self.name)
        seq = self.sequence
        wrapped = "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))
        return f">{hdr}\n{wrapped}\n"


@dataclass
class TraceAsset(AssetBase):
    sequence: str | None = None
    qualities: list[int] | None = None  # per-base PHRED, if present
    base_positions: list[int] | None = None  # PLOC2 (peak indices), if present
    channels: dict[BaseChar, list[int] | None] = field(
        default_factory=lambda: {"A": None, "C": None, "G": None, "T": None}
    )
    meta: dict[str, Any] = field(default_factory=dict)  # subset of ABIF tags, etc.

    def __post_init__(self) -> None:
        self.kind = AssetKind.TRACE

    def to_fasta(self, header: str | None = None) -> str | None:
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
    asset_ids: list[str] = field(default_factory=list)
    description: str = ""
    tags: list[str] = field(default_factory=list)
    primary_asset_id: str | None = None
    sequence_override: str | None = None
    feature_overrides: list[dict[str, Any]] | None = None
    created_at: datetime = field(default_factory=lambda: datetime.now(UTC))

    def effective_sequence(self, assets: dict[str, AssetBase]) -> str | None:
        if self.sequence_override:
            return self.sequence_override
        cand_id = self.primary_asset_id or next(
            (aid for aid in self.asset_ids if aid in assets), None
        )
        if not cand_id:
            return None
        a = assets.get(cand_id)
        if isinstance(a, SequenceAsset):
            return a.sequence
        if isinstance(a, TraceAsset):
            return a.sequence or None
        return None

    def effective_features(self, assets: dict[str, AssetBase]) -> list[dict[str, Any]]:
        if self.feature_overrides is not None:
            return self.feature_overrides
        cand_id = self.primary_asset_id or next(
            (aid for aid in self.asset_ids if aid in assets), None
        )
        a = assets.get(cand_id) if cand_id else None
        if isinstance(a, SequenceAsset) and a.features:
            return a.features
        return []

    def to_fasta(self, assets: dict[str, AssetBase]) -> str | None:
        seq = self.effective_sequence(assets)
        if not seq:
            return None
        hdr = self.name
        wrapped = "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))
        return f">{hdr}\n{wrapped}\n"
