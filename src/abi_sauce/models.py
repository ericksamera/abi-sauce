from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any
from datetime import datetime

class AssetKind(str, Enum):
    SEQUENCE = "sequence"      # e.g., FASTA/GenBank/ApE
    TRACE = "trace"            # e.g., AB1 (Sanger chromatogram)
    

@dataclass
class AssetBase:
    id: str
    name: str
    ext: str
    size: int
    checksum_md5: str
    created_at: datetime = field(default_factory=datetime.utcnow, init=False)  # <-- change
    kind: AssetKind = field(init=False)


@dataclass
class SequenceAsset(AssetBase):
    sequence: str               # plain sequence (A/C/G/T, possibly N)
    description: str = ""
    length: int = 0
    meta: Dict[str, Any] = field(default_factory=dict)   # e.g., source, locus, etc.
    features: List[Dict[str, Any]] = field(default_factory=list)  # simplified feature dicts

    def __post_init__(self):
        self.kind = AssetKind.SEQUENCE
        if not self.length:
            self.length = len(self.sequence)

    def to_fasta(self, header: Optional[str] = None) -> str:
        hdr = header or (self.description or self.name)
        seq = self.sequence
        # wrap 70-chars per FASTA convention
        wrapped = "\n".join(seq[i:i+70] for i in range(0, len(seq), 70))
        return f">{hdr}\n{wrapped}\n"

@dataclass
class TraceAsset(AssetBase):
    sequence: Optional[str] = None
    qualities: Optional[List[int]] = None    # per-base PHRED, if present
    base_positions: Optional[List[int]] = None  # peak indices (PLOC2), if present
    channels: Dict[str, Optional[List[int]]] = field(default_factory=lambda: {"A": None, "C": None, "G": None, "T": None})
    meta: Dict[str, Any] = field(default_factory=dict)   # e.g., instrument, run date, ABI tags subset

    def __post_init__(self):
        self.kind = AssetKind.TRACE

    def to_fasta(self, header: Optional[str] = None) -> Optional[str]:
        if not self.sequence:
            return None
        hdr = header or self.name
        wrapped = "\n".join(self.sequence[i:i+70] for i in range(0, len(self.sequence), 70))
        return f">{hdr}\n{wrapped}\n"