from __future__ import annotations

from abi_sauce import assembly as assembly_module
from abi_sauce import assembly_pairwise


def test_assembly_module_reexports_pairwise_entrypoints() -> None:
    assert (
        assembly_module.build_assembly_aligner
        is assembly_pairwise.build_assembly_aligner
    )
    assert (
        assembly_module.assemble_trimmed_pair is assembly_pairwise.assemble_trimmed_pair
    )
