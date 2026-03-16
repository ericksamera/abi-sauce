"""Compatibility facade for legacy assembly-service imports."""

from __future__ import annotations

from abi_sauce.services.assembly_compute import (
    AssemblyComputationStatus,
    ComputedAssembly,
    compute_saved_assemblies,
    compute_saved_assembly,
    compute_saved_multi_assembly,
)
from abi_sauce.services.assembly_export import (
    AssemblyDownloadArtifact,
    AssemblyExportFormat,
    AssemblyExportSelection,
    accepted_consensus_records,
    prepare_assembly_download,
    select_assembly_export,
)

__all__ = [
    "AssemblyComputationStatus",
    "AssemblyDownloadArtifact",
    "AssemblyExportFormat",
    "AssemblyExportSelection",
    "ComputedAssembly",
    "accepted_consensus_records",
    "compute_saved_assemblies",
    "compute_saved_assembly",
    "compute_saved_multi_assembly",
    "prepare_assembly_download",
    "select_assembly_export",
]
