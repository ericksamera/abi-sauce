from __future__ import annotations

from abi_sauce import assembly as assembly_module
from abi_sauce import assembly_multi


def test_assembly_module_reexports_multi_entrypoint() -> None:
    assert (
        assembly_module.assemble_trimmed_multi is assembly_multi.assemble_trimmed_multi
    )
