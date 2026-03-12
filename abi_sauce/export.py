from __future__ import annotations

from abi_sauce.exceptions import ExportError
from abi_sauce.models import SequenceRecord


def to_fasta(record: SequenceRecord, *, line_width: int = 80) -> str:
    """Serialize a sequence record as FASTA text."""
    if line_width <= 0:
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


def _header_for_record(record: SequenceRecord) -> str:
    """Build a stable single-line export header for a sequence record."""
    header = record.name or record.record_id or "sequence"
    return " ".join(header.split()) or "sequence"


def _wrap_sequence(sequence: str, line_width: int) -> str:
    """Wrap a sequence string to a fixed width for FASTA output."""
    if not sequence:
        return ""

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
