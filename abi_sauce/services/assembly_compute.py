from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from abi_sauce.assembly_exports import (
    consensus_record_from_multi_result,
    consensus_record_from_result,
)
from abi_sauce.assembly_multi import assemble_trimmed_multi
from abi_sauce.assembly_pairwise import assemble_trimmed_pair
from abi_sauce.assembly_state import AssemblyDefinition
from abi_sauce.assembly_types import AssemblyComputationResult
from abi_sauce.models import SequenceRecord
from abi_sauce.services.batch_trim import PreparedBatch

AssemblyComputationStatus = Literal[
    "ok",
    "invalid_definition",
    "not_implemented",
    "rejected",
    "error",
]


@dataclass(frozen=True, slots=True)
class ComputedAssembly:
    """One saved assembly definition resolved against the current prepared batch."""

    definition: AssemblyDefinition
    result: AssemblyComputationResult | None
    status: AssemblyComputationStatus
    status_reason: str | None = None
    consensus_record: SequenceRecord | None = None

    @property
    def is_accepted(self) -> bool:
        return self.status == "ok"

    @property
    def is_exportable(self) -> bool:
        return self.consensus_record is not None


def compute_saved_assembly(
    prepared_batch: PreparedBatch,
    definition: AssemblyDefinition,
) -> ComputedAssembly:
    """Resolve one saved assembly definition against the current prepared batch."""
    definition_reasons = _definition_ineligible_reasons(prepared_batch, definition)
    if definition_reasons:
        return ComputedAssembly(
            definition=definition,
            result=None,
            status="invalid_definition",
            status_reason="; ".join(definition_reasons),
            consensus_record=None,
        )

    if definition.engine_kind == "pairwise":
        return _compute_saved_pairwise_assembly(prepared_batch, definition)
    if definition.engine_kind == "multi":
        return compute_saved_multi_assembly(prepared_batch, definition)

    return ComputedAssembly(
        definition=definition,
        result=None,
        status="invalid_definition",
        status_reason=f"unsupported assembly engine: {definition.engine_kind}",
        consensus_record=None,
    )


def _compute_saved_pairwise_assembly(
    prepared_batch: PreparedBatch,
    definition: AssemblyDefinition,
) -> ComputedAssembly:
    pairwise_reasons = _pairwise_definition_ineligible_reasons(definition)
    if pairwise_reasons:
        return ComputedAssembly(
            definition=definition,
            result=None,
            status="invalid_definition",
            status_reason="; ".join(pairwise_reasons),
            consensus_record=None,
        )

    left_source_filename, right_source_filename = definition.source_filenames
    left_raw_record = prepared_batch.parsed_records[left_source_filename]
    right_raw_record = prepared_batch.parsed_records[right_source_filename]
    left_trim_result = prepared_batch.trim_results[left_source_filename]
    right_trim_result = prepared_batch.trim_results[right_source_filename]

    try:
        result = assemble_trimmed_pair(
            left_source_filename=left_source_filename,
            left_raw_record=left_raw_record,
            left_trim_result=left_trim_result,
            right_source_filename=right_source_filename,
            right_raw_record=right_raw_record,
            right_trim_result=right_trim_result,
            config=definition.config,
        )
    except Exception as exc:  # pragma: no cover - defensive boundary
        return ComputedAssembly(
            definition=definition,
            result=None,
            status="error",
            status_reason=str(exc),
            consensus_record=None,
        )

    consensus_record = (
        consensus_record_from_result(result, name=definition.name)
        if result.consensus_sequence
        else None
    )
    if result.accepted:
        return ComputedAssembly(
            definition=definition,
            result=result,
            status="ok",
            status_reason=None,
            consensus_record=consensus_record,
        )

    return ComputedAssembly(
        definition=definition,
        result=result,
        status="rejected",
        status_reason=result.rejection_reason or "assembly rejected",
        consensus_record=consensus_record,
    )


def compute_saved_multi_assembly(
    prepared_batch: PreparedBatch,
    definition: AssemblyDefinition,
) -> ComputedAssembly:
    """Resolve one saved multi-read assembly definition against the batch."""
    multi_reasons = _multi_definition_ineligible_reasons(definition)
    if multi_reasons:
        return ComputedAssembly(
            definition=definition,
            result=None,
            status="invalid_definition",
            status_reason="; ".join(multi_reasons),
            consensus_record=None,
        )

    try:
        result = assemble_trimmed_multi(
            source_filenames=definition.source_filenames,
            raw_records_by_source_filename=prepared_batch.parsed_records,
            trim_results_by_source_filename=prepared_batch.trim_results,
            config=definition.config,
        )
    except Exception as exc:  # pragma: no cover - defensive boundary
        return ComputedAssembly(
            definition=definition,
            result=None,
            status="error",
            status_reason=str(exc),
            consensus_record=None,
        )

    consensus_record = (
        consensus_record_from_multi_result(result, name=definition.name)
        if result.consensus_sequence
        else None
    )
    if result.accepted:
        return ComputedAssembly(
            definition=definition,
            result=result,
            status="ok",
            status_reason=None,
            consensus_record=consensus_record,
        )

    return ComputedAssembly(
        definition=definition,
        result=result,
        status="rejected",
        status_reason=result.rejection_reason or "assembly rejected",
        consensus_record=consensus_record,
    )


def compute_saved_assemblies(
    prepared_batch: PreparedBatch,
    definitions: Iterable[AssemblyDefinition],
) -> dict[str, ComputedAssembly]:
    """Resolve all saved assembly definitions against the current prepared batch."""
    return {
        definition.assembly_id: compute_saved_assembly(prepared_batch, definition)
        for definition in definitions
    }


def _definition_ineligible_reasons(
    prepared_batch: PreparedBatch,
    definition: AssemblyDefinition,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if definition.engine_kind not in {"pairwise", "multi"}:
        reasons.append(f"unsupported assembly engine: {definition.engine_kind}")
    if not definition.source_filenames:
        reasons.append("assembly requires at least 1 read")
    if len(set(definition.source_filenames)) != len(definition.source_filenames):
        reasons.append("assembly members must be distinct reads")

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

    return tuple(reasons)


def _pairwise_definition_ineligible_reasons(
    definition: AssemblyDefinition,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if len(definition.source_filenames) != 2:
        reasons.append("pairwise assembly currently requires exactly 2 reads")

    return tuple(reasons)


def _multi_definition_ineligible_reasons(
    definition: AssemblyDefinition,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if len(definition.source_filenames) < 2:
        reasons.append("multi assembly currently requires at least 2 reads")

    return tuple(reasons)


__all__ = [
    "AssemblyComputationStatus",
    "ComputedAssembly",
    "compute_saved_assemblies",
    "compute_saved_assembly",
    "compute_saved_multi_assembly",
]
