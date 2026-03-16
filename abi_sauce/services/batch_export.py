from __future__ import annotations

from dataclasses import dataclass

from abi_sauce.batch import ExportFormat
from abi_sauce.export import to_fasta_batch, to_fastq_batch, to_zip_batch
from abi_sauce.models import SequenceRecord
from abi_sauce.services.batch_trim import PreparedBatch


@dataclass(frozen=True, slots=True)
class BatchExportSelection:
    """Resolved export subset for one export choice."""

    export_format: ExportFormat
    require_min_length: bool
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]


@dataclass(frozen=True, slots=True)
class BatchDownloadArtifact:
    """Resolved export subset plus serialized batch download payload."""

    export_format: ExportFormat
    require_min_length: bool
    concatenate_batch: bool
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]
    data: str | bytes = ""
    filename: str = ""
    mime: str = ""

    @property
    def is_downloadable(self) -> bool:
        """Return whether a concrete download payload is available."""
        return bool(self.filename)


def select_batch_export(
    prepared_batch: PreparedBatch,
    *,
    export_format: ExportFormat,
    require_min_length: bool = False,
) -> BatchExportSelection:
    """Resolve the currently eligible batch export subset."""
    return BatchExportSelection(
        export_format=export_format,
        require_min_length=require_min_length,
        eligible_records=prepared_batch.batch_export_policy.eligible_records(
            export_format=export_format,
            require_min_length=require_min_length,
        ),
        ineligible_reasons=prepared_batch.batch_export_policy.ineligible_reasons_by_filename(
            export_format=export_format,
            require_min_length=require_min_length,
        ),
    )


def prepare_batch_download(
    prepared_batch: PreparedBatch,
    *,
    export_format: ExportFormat,
    concatenate_batch: bool,
    filename_stem: str,
    require_min_length: bool = False,
    fasta_line_width: int | None = 80,
) -> BatchDownloadArtifact:
    """Build one batch download artifact for the current export selection."""
    export_selection = select_batch_export(
        prepared_batch,
        export_format=export_format,
        require_min_length=require_min_length,
    )

    if not export_selection.eligible_records:
        return BatchDownloadArtifact(
            export_format=export_format,
            require_min_length=require_min_length,
            concatenate_batch=concatenate_batch,
            eligible_records=export_selection.eligible_records,
            ineligible_reasons=export_selection.ineligible_reasons,
        )

    export_data, export_filename, export_mime = _serialize_batch_download(
        records=export_selection.eligible_records,
        export_format=export_format,
        concatenate_batch=concatenate_batch,
        filename_stem=filename_stem,
        fasta_line_width=fasta_line_width,
    )

    return BatchDownloadArtifact(
        export_format=export_format,
        require_min_length=require_min_length,
        concatenate_batch=concatenate_batch,
        eligible_records=export_selection.eligible_records,
        ineligible_reasons=export_selection.ineligible_reasons,
        data=export_data,
        filename=export_filename,
        mime=export_mime,
    )


def _serialize_batch_download(
    *,
    records: tuple[SequenceRecord, ...],
    export_format: ExportFormat,
    concatenate_batch: bool,
    filename_stem: str,
    fasta_line_width: int | None = 80,
) -> tuple[str | bytes, str, str]:
    """Serialize the current eligible batch export selection."""
    normalized_filename_stem = filename_stem.strip() or "abi-sauce-batch"

    if concatenate_batch and export_format == "fasta":
        return (
            to_fasta_batch(records, line_width=fasta_line_width),
            f"{normalized_filename_stem}.fasta",
            "text/plain",
        )

    if concatenate_batch and export_format == "fastq":
        return (
            to_fastq_batch(records),
            f"{normalized_filename_stem}.fastq",
            "text/plain",
        )

    return (
        to_zip_batch(
            records,
            export_format=export_format,
            line_width=fasta_line_width,
        ),
        f"{normalized_filename_stem}.zip",
        "application/zip",
    )


__all__ = [
    "BatchDownloadArtifact",
    "BatchExportSelection",
    "prepare_batch_download",
    "select_batch_export",
]
