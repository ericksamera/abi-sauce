from __future__ import annotations
from typing import List, Optional
from abi_sauce.models import TraceAsset, SequenceAsset

def consensus_from_traces(traces: List[TraceAsset], name: Optional[str] = None) -> SequenceAsset:
    """
    Placeholder for future electropherogram alignment + consensus calling.

    Expected behavior later:
    - basecall normalization & trimming
    - pairwise/multiple alignment of called bases (plus peak-aware tie-breaks)
    - quality-aware consensus generation -> SequenceAsset
    """
    raise NotImplementedError("Consensus calling not implemented yet.")
