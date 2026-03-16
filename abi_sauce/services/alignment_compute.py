from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from abi_sauce.alignment_state import AlignmentDefinition
from abi_sauce.assembly_state import AssemblyDefinition as SavedAssemblyDefinition
from abi_sauce.services.assembly_compute import (
    ComputedAssembly,
    compute_saved_assembly,
)
from abi_sauce.services.batch_trim import PreparedBatch
from abi_sauce.services.reference_alignment import (
    ComputedReferenceAlignment,
    ComputedReferenceMultiAlignment,
    compute_reference_alignment,
    compute_reference_multi_alignment,
)

AlignmentComputationStatus = Literal[
    "ok",
    "invalid_definition",
    "not_implemented",
    "rejected",
    "error",
]


@dataclass(frozen=True, slots=True)
class ComputedAlignment:
    """One saved alignment definition resolved against the current prepared batch."""

    definition: AlignmentDefinition
    status: AlignmentComputationStatus
    status_reason: str | None = None
    assembly: ComputedAssembly | None = None
    reference_alignment: ComputedReferenceAlignment | None = None
    reference_multi_alignment: ComputedReferenceMultiAlignment | None = None

    @property
    def is_ready(self) -> bool:
        return (
            self.assembly is not None
            or self.reference_alignment is not None
            or self.reference_multi_alignment is not None
        )

    @property
    def is_accepted(self) -> bool:
        return self.status == "ok"


def compute_saved_alignment(
    prepared_batch: PreparedBatch,
    definition: AlignmentDefinition,
) -> ComputedAlignment:
    """Resolve one saved alignment definition against the current prepared batch."""
    definition_reasons = _definition_ineligible_reasons(prepared_batch, definition)
    if definition_reasons:
        return ComputedAlignment(
            definition=definition,
            status="invalid_definition",
            status_reason="; ".join(definition_reasons),
        )

    if definition.engine_kind == "pairwise":
        return _compute_saved_pairwise_alignment(prepared_batch, definition)
    if definition.engine_kind == "reference_single":
        return _compute_saved_reference_alignment(prepared_batch, definition)
    if definition.engine_kind == "reference_multi":
        return _compute_saved_reference_multi_alignment(prepared_batch, definition)

    return ComputedAlignment(
        definition=definition,
        status="invalid_definition",
        status_reason=f"unsupported alignment engine: {definition.engine_kind}",
    )


def compute_saved_alignments(
    prepared_batch: PreparedBatch,
    definitions: Iterable[AlignmentDefinition],
) -> dict[str, ComputedAlignment]:
    """Resolve all saved alignment definitions against the current prepared batch."""
    return {
        definition.alignment_id: compute_saved_alignment(prepared_batch, definition)
        for definition in definitions
    }


def _compute_saved_pairwise_alignment(
    prepared_batch: PreparedBatch,
    definition: AlignmentDefinition,
) -> ComputedAlignment:
    computed_assembly = compute_saved_assembly(
        prepared_batch,
        SavedAssemblyDefinition(
            assembly_id=definition.alignment_id,
            name=definition.name,
            source_filenames=definition.source_filenames,
            config=definition.assembly_config,
            engine_kind="pairwise",
        ),
    )
    return ComputedAlignment(
        definition=definition,
        status=computed_assembly.status,
        status_reason=computed_assembly.status_reason,
        assembly=computed_assembly,
    )


def _compute_saved_reference_alignment(
    prepared_batch: PreparedBatch,
    definition: AlignmentDefinition,
) -> ComputedAlignment:
    try:
        computed_reference_alignment = compute_reference_alignment(
            prepared_batch,
            source_filename=definition.source_filenames[0],
            reference_text=definition.reference_text or "",
            strand_policy=definition.strand_policy,
        )
    except ValueError as exc:
        return ComputedAlignment(
            definition=definition,
            status="error",
            status_reason=str(exc),
        )
    except Exception as exc:  # pragma: no cover - defensive boundary
        return ComputedAlignment(
            definition=definition,
            status="error",
            status_reason=str(exc),
        )

    return ComputedAlignment(
        definition=definition,
        status="ok",
        status_reason=None,
        reference_alignment=computed_reference_alignment,
    )


def _compute_saved_reference_multi_alignment(
    prepared_batch: PreparedBatch,
    definition: AlignmentDefinition,
) -> ComputedAlignment:
    try:
        computed_reference_multi_alignment = compute_reference_multi_alignment(
            prepared_batch,
            source_filenames=definition.source_filenames,
            reference_text=definition.reference_text or "",
            reference_name=definition.reference_name,
            strand_policy=definition.strand_policy,
            config=definition.assembly_config,
        )
    except ValueError as exc:
        return ComputedAlignment(
            definition=definition,
            status="error",
            status_reason=str(exc),
        )
    except Exception as exc:  # pragma: no cover - defensive boundary
        return ComputedAlignment(
            definition=definition,
            status="error",
            status_reason=str(exc),
        )

    if computed_reference_multi_alignment.result.accepted:
        return ComputedAlignment(
            definition=definition,
            status="ok",
            status_reason=None,
            reference_multi_alignment=computed_reference_multi_alignment,
        )

    return ComputedAlignment(
        definition=definition,
        status="rejected",
        status_reason=(
            computed_reference_multi_alignment.result.rejection_reason
            or "reference-guided multi-read alignment rejected"
        ),
        reference_multi_alignment=computed_reference_multi_alignment,
    )


def _definition_ineligible_reasons(
    prepared_batch: PreparedBatch,
    definition: AlignmentDefinition,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if definition.engine_kind not in {
        "pairwise",
        "reference_single",
        "reference_multi",
    }:
        reasons.append(f"unsupported alignment engine: {definition.engine_kind}")
    if not definition.source_filenames:
        reasons.append("alignment requires at least 1 read")
    if len(set(definition.source_filenames)) != len(definition.source_filenames):
        reasons.append("alignment members must be distinct reads")

    missing_filenames = [
        source_filename
        for source_filename in definition.source_filenames
        if source_filename not in prepared_batch.parsed_records
        or source_filename not in prepared_batch.trim_results
    ]
    if missing_filenames:
        reasons.append(
            "missing active batch members: " + ", ".join(sorted(missing_filenames))
        )

    if definition.engine_kind == "pairwise" and len(definition.source_filenames) != 2:
        reasons.append("pairwise alignment currently requires exactly 2 reads")
    if (
        definition.engine_kind == "reference_single"
        and len(definition.source_filenames) != 1
    ):
        reasons.append("single-reference alignment currently requires exactly 1 read")
    if (
        definition.engine_kind == "reference_multi"
        and len(definition.source_filenames) < 2
    ):
        reasons.append(
            "reference-guided multi-read alignment currently requires at least 2 reads"
        )
    if definition.engine_kind in {"reference_single", "reference_multi"} and (
        not isinstance(definition.reference_text, str)
        or not definition.reference_text.strip()
    ):
        reasons.append("reference-guided alignment requires reference text")

    return tuple(reasons)


__all__ = [
    "AlignmentComputationStatus",
    "ComputedAlignment",
    "compute_saved_alignment",
    "compute_saved_alignments",
]
