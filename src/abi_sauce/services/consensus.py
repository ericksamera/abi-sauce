from __future__ import annotations

from abi_sauce.models import SequenceAsset, TraceAsset


def consensus_from_traces(
    traces: list[TraceAsset], name: str | None = None
) -> SequenceAsset:
    """
    Placeholder for future electropherogram alignment + consensus calling.

    Expected behavior later:
    - basecall normalization & trimming
    - pairwise/multiple alignment of called bases (plus peak-aware tie-breaks)
    - quality-aware consensus generation -> SequenceAsset
    """
    raise NotImplementedError("Consensus calling not implemented yet.")
