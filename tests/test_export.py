from __future__ import annotations

import pytest

from abi_sauce.export import ExportError, to_fasta, to_fastq
from abi_sauce.models import SequenceRecord
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    sequence: str = "ACGTACGT",
    qualities: list[int] | None = None,
    name: str = "trace_001",
) -> SequenceRecord:
    return SequenceRecord(
        record_id="trace_001",
        name=name,
        description="test record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
    )


def test_to_fasta_serializes_header_and_sequence() -> None:
    record = make_record(sequence="ACGT")
    result = to_fasta(record)

    assert result == ">trace_001\nACGT\n"


def test_to_fasta_wraps_sequence() -> None:
    record = make_record(sequence="ACGTACGTAC")
    result = to_fasta(record, line_width=4)

    assert result == ">trace_001\nACGT\nACGT\nAC\n"


def test_to_fasta_rejects_invalid_line_width() -> None:
    record = make_record()

    with pytest.raises(ValueError, match="line_width must be > 0"):
        to_fasta(record, line_width=0)


def test_to_fastq_serializes_four_line_record() -> None:
    record = make_record(
        sequence="ACGT",
        qualities=[40, 41, 42, 43],
    )
    result = to_fastq(record)

    assert result == "@trace_001\nACGT\n+\nIJKL\n"


def test_to_fastq_requires_qualities() -> None:
    record = make_record(sequence="ACGT", qualities=None)

    with pytest.raises(ExportError, match="FASTQ export requires per-base qualities."):
        to_fastq(record)


def test_to_fastq_requires_matching_sequence_and_quality_lengths() -> None:
    record = make_record(
        sequence="ACGT",
        qualities=[40, 41, 42],
    )

    with pytest.raises(
        ExportError,
        match="FASTQ export requires sequence and qualities to have the same length.",
    ):
        to_fastq(record)


def test_to_fastq_rejects_out_of_range_phred_scores() -> None:
    record = make_record(
        sequence="ACGT",
        qualities=[40, 41, 42, 120],
    )

    with pytest.raises(
        ExportError,
        match=r"FASTQ export requires PHRED scores between 0 and 93, got: 120",
    ):
        to_fastq(record)


def test_trimmed_record_exports_trimmed_fasta_sequence() -> None:
    record = make_record(sequence="AACCGGTT")
    trim_result = trim_sequence_record(
        record,
        TrimConfig(left_trim=2, right_trim=2),
    )

    result = to_fasta(trim_result.record)
    assert result == ">trace_001\nCCGG\n"


def test_trimmed_record_exports_trimmed_fastq_sequence_and_qualities() -> None:
    record = make_record(
        sequence="AACCGGTT",
        qualities=[30, 31, 32, 33, 34, 35, 36, 37],
    )
    trim_result = trim_sequence_record(
        record,
        TrimConfig(left_trim=2, right_trim=2),
    )

    result = to_fastq(trim_result.record)
    assert result == "@trace_001\nCCGG\n+\nABCD\n"


def test_export_header_collapses_whitespace() -> None:
    record = make_record(name="trace   001\nsample", sequence="ACGT")
    result = to_fasta(record)

    assert result == ">trace 001 sample\nACGT\n"
