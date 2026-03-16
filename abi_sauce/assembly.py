"""Compatibility facade for legacy assembly imports."""

from __future__ import annotations

from abi_sauce.assembly_exports import (
    consensus_record_from_multi_result,
    consensus_record_from_result,
    format_assembly_alignment_fasta,
)
from abi_sauce.assembly_multi import assemble_trimmed_multi
from abi_sauce.assembly_pairwise import assemble_trimmed_pair, build_assembly_aligner
from abi_sauce.assembly_presenters import (
    assembly_conflicts_to_rows,
    format_assembly_block,
)
from abi_sauce.assembly_types import (
    AssemblyComputationResult,
    AssemblyConfig,
    AssemblyResult,
    MultiAssemblyResult,
)

__all__ = [
    "AssemblyComputationResult",
    "AssemblyConfig",
    "AssemblyResult",
    "MultiAssemblyResult",
    "build_assembly_aligner",
    "assemble_trimmed_pair",
    "assemble_trimmed_multi",
    "assembly_conflicts_to_rows",
    "format_assembly_block",
    "format_assembly_alignment_fasta",
    "consensus_record_from_result",
    "consensus_record_from_multi_result",
]
