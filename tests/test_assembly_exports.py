from __future__ import annotations

from abi_sauce import assembly as assembly_module
from abi_sauce import assembly_exports


def test_assembly_module_reexports_export_helpers() -> None:
    assert (
        assembly_module.format_assembly_alignment_fasta
        is assembly_exports.format_assembly_alignment_fasta
    )
    assert (
        assembly_module.consensus_record_from_result
        is assembly_exports.consensus_record_from_result
    )
    assert (
        assembly_module.consensus_record_from_multi_result
        is assembly_exports.consensus_record_from_multi_result
    )
