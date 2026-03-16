from __future__ import annotations

from abi_sauce import reference_alignment as reference_alignment_module
from abi_sauce import reference_alignment_presenters
from abi_sauce import reference_alignment_types


def test_reference_alignment_module_reexports_types_and_presenters() -> None:
    assert (
        reference_alignment_module.AlignmentEvent
        is reference_alignment_types.AlignmentEvent
    )
    assert (
        reference_alignment_module.AlignmentResult
        is reference_alignment_types.AlignmentResult
    )
    assert (
        reference_alignment_module.ReferenceAlignmentColumn
        is reference_alignment_types.ReferenceAlignmentColumn
    )
    assert (
        reference_alignment_module.alignment_events_to_rows
        is reference_alignment_presenters.alignment_events_to_rows
    )
    assert (
        reference_alignment_module.format_alignment_block
        is reference_alignment_presenters.format_alignment_block
    )
