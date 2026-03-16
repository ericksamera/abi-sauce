from __future__ import annotations

from abi_sauce import assembly as assembly_module
from abi_sauce import assembly_presenters


def test_assembly_module_reexports_presenter_helpers() -> None:
    assert (
        assembly_module.assembly_conflicts_to_rows
        is assembly_presenters.assembly_conflicts_to_rows
    )
    assert (
        assembly_module.format_assembly_block
        is assembly_presenters.format_assembly_block
    )
