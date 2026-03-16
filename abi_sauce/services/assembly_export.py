from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from io import BytesIO
import json
from typing import Literal
import zipfile

from abi_sauce.assembly_presenters import format_assembly_block
from abi_sauce.assembly_types import AssemblyResult, MultiAssemblyResult
from abi_sauce.exceptions import ExportError
from abi_sauce.export import to_fasta, to_fasta_batch, to_zip_batch
from abi_sauce.models import SequenceRecord
from abi_sauce.services.assembly_compute import ComputedAssembly

AssemblyExportFormat = Literal["fasta"]


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


def _computed_assembly_ineligible_reasons(
    computed_assembly: ComputedAssembly,
    *,
    require_accepted: bool,
) -> tuple[str, ...]:
    reasons: list[str] = []
    blocked_by_status = False

    if computed_assembly.status == "invalid_definition":
        reasons.append(
            computed_assembly.status_reason or "assembly definition is invalid"
        )
        blocked_by_status = True
    elif computed_assembly.status == "not_implemented":
        reasons.append(
            computed_assembly.status_reason or "assembly engine is not implemented yet"
        )
        blocked_by_status = True
    elif computed_assembly.status == "error":
        reasons.append(computed_assembly.status_reason or "assembly failed to compute")
        blocked_by_status = True
    elif require_accepted and computed_assembly.status == "rejected":
        reasons.append(computed_assembly.status_reason or "assembly was rejected")

    if computed_assembly.consensus_record is None and not blocked_by_status:
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
    manifest_entry: dict[str, object] = {
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
        "rejection_reason": None if result is None else result.rejection_reason,
    }
    if isinstance(result, AssemblyResult):
        manifest_entry.update(
            {
                "chosen_right_orientation": result.chosen_right_orientation,
                "overlap_length": result.overlap_length,
                "percent_identity": result.percent_identity,
                "conflict_count": result.conflict_count,
                "alignment_preview": format_assembly_block(result),
            }
        )
    elif isinstance(result, MultiAssemblyResult):
        manifest_entry.update(
            {
                "seed_member_index": result.seed_member_index,
                "included_member_count": result.included_member_count,
                "excluded_member_count": result.excluded_member_count,
                "ambiguous_column_count": result.ambiguous_column_count,
                "member_summaries": [
                    {
                        "member_index": member.member_index,
                        "source_filename": member.source_filename,
                        "display_name": member.display_name,
                        "is_seed": member.is_seed,
                        "included": member.included,
                        "chosen_orientation": member.chosen_orientation,
                        "inclusion_reason": member.inclusion_reason,
                        "overlap_length": member.overlap_length,
                        "percent_identity": member.percent_identity,
                    }
                    for member in result.members
                ],
                "gapped_consensus": result.gapped_consensus,
            }
        )
    return manifest_entry


def _safe_filename(value: str) -> str:
    collapsed = "_".join(value.split())
    return "".join(
        character
        for character in collapsed
        if character.isalnum() or character in {"-", "_", "."}
    ).strip("._-")


__all__ = [
    "AssemblyDownloadArtifact",
    "AssemblyExportFormat",
    "AssemblyExportSelection",
    "accepted_consensus_records",
    "prepare_assembly_download",
    "select_assembly_export",
]
