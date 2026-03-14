from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from io import BytesIO
import json
from typing import Literal
import zipfile

from abi_sauce.assembly import (
    AssemblyResult,
    assemble_trimmed_pair,
    consensus_record_from_result,
    format_assembly_block,
)
from abi_sauce.assembly_state import AssemblyDefinition
from abi_sauce.exceptions import ExportError
from abi_sauce.export import to_fasta, to_fasta_batch, to_zip_batch
from abi_sauce.models import SequenceRecord
from abi_sauce.services.batch import PreparedBatch

AssemblyComputationStatus = Literal["ok", "invalid_definition", "rejected", "error"]
AssemblyExportFormat = Literal["fasta"]


@dataclass(frozen=True, slots=True)
class ComputedAssembly:
    """One saved assembly definition resolved against the current prepared batch."""

    definition: AssemblyDefinition
    result: AssemblyResult | None
    status: AssemblyComputationStatus
    status_reason: str | None = None
    consensus_record: SequenceRecord | None = None

    @property
    def is_accepted(self) -> bool:
        return self.status == "ok"

    @property
    def is_exportable(self) -> bool:
        return self.consensus_record is not None


@dataclass(frozen=True, slots=True)
class AssemblyExportSelection:
    """Resolved consensus export subset for saved assemblies."""

    require_accepted: bool
    eligible_assemblies: tuple[ComputedAssembly, ...]
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]


@dataclass(frozen=True, slots=True)
class AssemblyDownloadArtifact:
    """Resolved saved-assembly export plus serialized download payload."""

    export_format: AssemblyExportFormat
    require_accepted: bool
    concatenate_batch: bool
    include_manifest: bool
    eligible_assemblies: tuple[ComputedAssembly, ...]
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]
    data: str | bytes = ""
    filename: str = ""
    mime: str = ""

    @property
    def is_downloadable(self) -> bool:
        return bool(self.filename)


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


def compute_saved_assemblies(
    prepared_batch: PreparedBatch,
    definitions: Iterable[AssemblyDefinition],
) -> dict[str, ComputedAssembly]:
    """Resolve all saved assembly definitions against the current prepared batch."""
    return {
        definition.assembly_id: compute_saved_assembly(prepared_batch, definition)
        for definition in definitions
    }


def accepted_consensus_records(
    computed_assemblies: Mapping[str, ComputedAssembly],
    *,
    selected_ids: Iterable[str] | None = None,
) -> tuple[SequenceRecord, ...]:
    """Return accepted exportable consensus records in display order."""
    export_selection = select_assembly_export(
        computed_assemblies,
        selected_ids=selected_ids,
        require_accepted=True,
    )
    return export_selection.eligible_records


def select_assembly_export(
    computed_assemblies: Mapping[str, ComputedAssembly],
    *,
    selected_ids: Iterable[str] | None = None,
    require_accepted: bool = True,
) -> AssemblyExportSelection:
    """Resolve which saved assembly consensus records are export-eligible."""
    selected_ids_tuple = (
        tuple(selected_ids) if selected_ids is not None else tuple(computed_assemblies)
    )
    selected_id_set = frozenset(selected_ids_tuple)

    eligible_assemblies: list[ComputedAssembly] = []
    ineligible_reasons: list[tuple[str, tuple[str, ...]]] = []

    for assembly_id, computed_assembly in computed_assemblies.items():
        if assembly_id not in selected_id_set:
            continue
        reasons = _computed_assembly_ineligible_reasons(
            computed_assembly,
            require_accepted=require_accepted,
        )
        if reasons:
            ineligible_reasons.append((_assembly_label(computed_assembly), reasons))
            continue
        eligible_assemblies.append(computed_assembly)

    for selected_id in selected_ids_tuple:
        if selected_id not in computed_assemblies:
            ineligible_reasons.append(
                (selected_id, ("assembly definition is no longer available",))
            )

    return AssemblyExportSelection(
        require_accepted=require_accepted,
        eligible_assemblies=tuple(eligible_assemblies),
        eligible_records=tuple(
            computed_assembly.consensus_record
            for computed_assembly in eligible_assemblies
            if computed_assembly.consensus_record is not None
        ),
        ineligible_reasons=tuple(ineligible_reasons),
    )


def prepare_assembly_download(
    computed_assemblies: Mapping[str, ComputedAssembly],
    *,
    selected_ids: Iterable[str] | None,
    concatenate_batch: bool,
    filename_stem: str,
    export_format: AssemblyExportFormat = "fasta",
    require_accepted: bool = True,
    include_manifest: bool = False,
    fasta_line_width: int | None = 80,
) -> AssemblyDownloadArtifact:
    """Build one bulk consensus download artifact from saved assemblies."""
    if export_format != "fasta":
        raise ExportError("Assembly export currently supports FASTA only.")

    export_selection = select_assembly_export(
        computed_assemblies,
        selected_ids=selected_ids,
        require_accepted=require_accepted,
    )
    if not export_selection.eligible_records:
        return AssemblyDownloadArtifact(
            export_format=export_format,
            require_accepted=require_accepted,
            concatenate_batch=concatenate_batch,
            include_manifest=include_manifest,
            eligible_assemblies=export_selection.eligible_assemblies,
            eligible_records=export_selection.eligible_records,
            ineligible_reasons=export_selection.ineligible_reasons,
        )

    normalized_filename_stem = filename_stem.strip() or "abi-sauce-assemblies"
    if concatenate_batch:
        return AssemblyDownloadArtifact(
            export_format=export_format,
            require_accepted=require_accepted,
            concatenate_batch=concatenate_batch,
            include_manifest=False,
            eligible_assemblies=export_selection.eligible_assemblies,
            eligible_records=export_selection.eligible_records,
            ineligible_reasons=export_selection.ineligible_reasons,
            data=to_fasta_batch(
                export_selection.eligible_records,
                line_width=fasta_line_width,
            ),
            filename=f"{normalized_filename_stem}.fasta",
            mime="text/plain",
        )

    zip_payload = (
        _to_zip_batch_with_manifest(
            export_selection.eligible_assemblies,
            line_width=fasta_line_width,
            ineligible_reasons=export_selection.ineligible_reasons,
        )
        if include_manifest
        else to_zip_batch(
            export_selection.eligible_records,
            export_format="fasta",
            line_width=fasta_line_width,
        )
    )
    return AssemblyDownloadArtifact(
        export_format=export_format,
        require_accepted=require_accepted,
        concatenate_batch=concatenate_batch,
        include_manifest=include_manifest,
        eligible_assemblies=export_selection.eligible_assemblies,
        eligible_records=export_selection.eligible_records,
        ineligible_reasons=export_selection.ineligible_reasons,
        data=zip_payload,
        filename=f"{normalized_filename_stem}.zip",
        mime="application/zip",
    )


def _definition_ineligible_reasons(
    prepared_batch: PreparedBatch,
    definition: AssemblyDefinition,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if definition.engine_kind != "pairwise":
        reasons.append(f"unsupported assembly engine: {definition.engine_kind}")
    if len(definition.source_filenames) != 2:
        reasons.append("pairwise assembly currently requires exactly 2 reads")
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


def _computed_assembly_ineligible_reasons(
    computed_assembly: ComputedAssembly,
    *,
    require_accepted: bool,
) -> tuple[str, ...]:
    reasons: list[str] = []

    if computed_assembly.status == "invalid_definition":
        reasons.append(
            computed_assembly.status_reason or "assembly definition is invalid"
        )
    elif computed_assembly.status == "error":
        reasons.append(computed_assembly.status_reason or "assembly failed to compute")
    elif require_accepted and computed_assembly.status == "rejected":
        reasons.append(computed_assembly.status_reason or "assembly was rejected")

    if computed_assembly.consensus_record is None:
        reasons.append("no consensus sequence is available")

    return tuple(reasons)


def _assembly_label(computed_assembly: ComputedAssembly) -> str:
    return computed_assembly.definition.name or computed_assembly.definition.assembly_id


def _to_zip_batch_with_manifest(
    eligible_assemblies: tuple[ComputedAssembly, ...],
    *,
    line_width: int | None,
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...],
) -> bytes:
    with BytesIO() as file_buffer:
        with zipfile.ZipFile(
            file_buffer,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_file:
            for index, computed_assembly in enumerate(eligible_assemblies, start=1):
                consensus_record = computed_assembly.consensus_record
                if consensus_record is None:
                    continue
                basename = _safe_filename(
                    consensus_record.name or consensus_record.record_id
                )
                if not basename:
                    basename = f"assembly_{index:03d}"
                zip_file.writestr(
                    f"{index:03d}_{basename}.fasta",
                    to_fasta(consensus_record, line_width=line_width),
                )
            zip_file.writestr(
                "manifest.json",
                json.dumps(
                    {
                        "exported_assemblies": [
                            _manifest_entry(computed_assembly)
                            for computed_assembly in eligible_assemblies
                        ],
                        "excluded_assemblies": [
                            {
                                "assembly_name": assembly_name,
                                "reasons": list(reasons),
                            }
                            for assembly_name, reasons in ineligible_reasons
                        ],
                    },
                    indent=2,
                    sort_keys=True,
                ),
            )
        return file_buffer.getvalue()


def _manifest_entry(computed_assembly: ComputedAssembly) -> dict[str, object]:
    result = computed_assembly.result
    definition = computed_assembly.definition
    return {
        "assembly_id": definition.assembly_id,
        "assembly_name": definition.name,
        "engine_kind": definition.engine_kind,
        "member_filenames": list(definition.source_filenames),
        "accepted": computed_assembly.is_accepted,
        "status": computed_assembly.status,
        "status_reason": computed_assembly.status_reason,
        "config": {
            "match_score": definition.config.match_score,
            "mismatch_score": definition.config.mismatch_score,
            "open_internal_gap_score": definition.config.open_internal_gap_score,
            "extend_internal_gap_score": definition.config.extend_internal_gap_score,
            "min_overlap_length": definition.config.min_overlap_length,
            "min_percent_identity": definition.config.min_percent_identity,
            "quality_margin": definition.config.quality_margin,
        },
        "consensus_record_name": (
            None
            if computed_assembly.consensus_record is None
            else computed_assembly.consensus_record.name
        ),
        "consensus_length": (
            None
            if computed_assembly.consensus_record is None
            else len(computed_assembly.consensus_record.sequence)
        ),
        "chosen_right_orientation": (
            None if result is None else result.chosen_right_orientation
        ),
        "overlap_length": None if result is None else result.overlap_length,
        "percent_identity": None if result is None else result.percent_identity,
        "conflict_count": None if result is None else result.conflict_count,
        "rejection_reason": None if result is None else result.rejection_reason,
        "alignment_preview": None if result is None else format_assembly_block(result),
    }


def _safe_filename(value: str) -> str:
    collapsed = "_".join(value.split())
    return "".join(
        character
        for character in collapsed
        if character.isalnum() or character in {"-", "_", "."}
    ).strip("._-")
