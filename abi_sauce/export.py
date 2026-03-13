from __future__ import annotations

from collections.abc import Iterable
from io import BytesIO
from typing import Literal
import zipfile

from abi_sauce.exceptions import ExportError
from abi_sauce.models import SequenceRecord


def to_fasta(record: SequenceRecord, *, line_width: int | None = 80) -> str:
    """Serialize a sequence record as FASTA text."""
    if line_width is not None and line_width <= 0:
        raise ValueError("line_width must be > 0")

    header = _header_for_record(record)
    wrapped_sequence = _wrap_sequence(record.sequence, line_width)

    return f">{header}\n{wrapped_sequence}\n"


def to_fastq(record: SequenceRecord) -> str:
    """Serialize a sequence record as FASTQ text using PHRED+33 encoding."""
    header = _header_for_record(record)
    qualities = record.qualities

    if qualities is None:
        raise ExportError("FASTQ export requires per-base qualities.")

    if len(qualities) != len(record.sequence):
        raise ExportError(
            "FASTQ export requires sequence and qualities to have the same length."
        )

    quality_string = "".join(_encode_phred_quality(score) for score in qualities)
    return f"@{header}\n{record.sequence}\n+\n{quality_string}\n"


def to_fasta_batch(
    records: Iterable[SequenceRecord],
    *,
    line_width: int | None = 80,
) -> str:
    """Serialize multiple sequence records into one FASTA string."""
    return "".join(to_fasta(record, line_width=line_width) for record in records)


def to_fastq_batch(records: Iterable[SequenceRecord]) -> str:
    """Serialize multiple sequence records into one FASTQ string."""
    return "".join(to_fastq(record) for record in records)


def to_zip_batch(
    records: Iterable[SequenceRecord],
    *,
    export_format: Literal["fasta", "fastq"],
    line_width: int | None = 80,
) -> bytes:
    """Serialize records into a ZIP of per-record FASTA/FASTQ files."""
    with BytesIO() as file_buffer:
        with zipfile.ZipFile(
            file_buffer,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_file:
            for index, record in enumerate(records, start=1):
                basename = _safe_filename(_header_for_record(record))
                if not basename:
                    basename = f"sequence_{index:03d}"

                zip_file.writestr(
                    f"{index:03d}_{basename}.{_extension_for_format(export_format)}",
                    _serialize_record(
                        record,
                        export_format=export_format,
                        line_width=line_width,
                    ),
                )

        return file_buffer.getvalue()


def _header_for_record(record: SequenceRecord) -> str:
    """Build a stable single-line export header for a sequence record."""
    header = record.name or record.record_id or "sequence"
    return " ".join(header.split()) or "sequence"


def _serialize_record(
    record: SequenceRecord,
    *,
    export_format: Literal["fasta", "fastq"],
    line_width: int | None = 80,
) -> str:
    """Serialize a single record in the requested format."""
    if export_format == "fasta":
        return to_fasta(record, line_width=line_width)
    if export_format == "fastq":
        return to_fastq(record)

    raise ExportError(f"Unsupported export format: {export_format}")


def _extension_for_format(export_format: Literal["fasta", "fastq"]) -> str:
    """Return the file extension for an export format."""
    return export_format


def _safe_filename(value: str) -> str:
    """Convert a record header into a filesystem-safe stem."""
    collapsed = "_".join(value.split())
    return "".join(
        char for char in collapsed if char.isalnum() or char in {"-", "_", "."}
    ).strip("._-")


def _wrap_sequence(sequence: str, line_width: int | None) -> str:
    """Wrap a sequence string to a fixed width for FASTA output."""
    if not sequence:
        return ""

    if line_width is None:
        return sequence

    return "\n".join(
        sequence[index : index + line_width]
        for index in range(0, len(sequence), line_width)
    )


def _encode_phred_quality(score: int) -> str:
    """Convert a PHRED score to a single FASTQ quality character."""
    if not 0 <= score <= 93:
        raise ExportError(
            f"FASTQ export requires PHRED scores between 0 and 93, got: {score}"
        )

    return chr(score + 33)
