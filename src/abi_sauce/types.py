# src/abi_sauce/types.py
from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, Protocol, TypedDict

# Canonical DNA base symbol (uppercase only)
BaseChar = Literal["A", "C", "G", "T"]

# Index of a called base in a read (0-based inclusive)
BaseIndex = int
# Inclusive sample-space window [start, end] in trace sample indices
SampleSpan = tuple[int, int]

# Numeric type accepted by many functions
Number = int | float

# Channel arrays are per-sample fluorescence intensities
Channel = Sequence(Number)
ChannelsMap = dict[BaseChar, Channel]

# Aggregated (per-base) values per channel
PerBase = Sequence(Number)
PerBaseMap = dict[BaseChar, PerBase]


class ABIFMeta(TypedDict, total=False):
    instrument: str
    run_date: str
    abif_order: str
    abif_keys: list[str]


class HasToFasta(Protocol):
    def to_fasta(self, header: str | None = None) -> str | None: ...


__all__ = [
    "BaseChar",
    "BaseIndex",
    "SampleSpan",
    "Number",
    "Channel",
    "ChannelsMap",
    "PerBase",
    "PerBaseMap",
    "ABIFMeta",
    "HasToFasta",
]
