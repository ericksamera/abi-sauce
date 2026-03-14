from __future__ import annotations

from abi_sauce.batch import build_batch_export_policy, build_batch_summary
from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData
from abi_sauce.trimming import TrimConfig, TrimResult, trim_sequence_record


def make_upload(filename: str) -> SequenceUpload:
    return SequenceUpload(filename=filename, content=b"fake")


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    trace_data: TraceData | None = None,
) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="test record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=trace_data,
        annotations={"source_filename": f"{name}.ab1"},
    )


def test_build_batch_summary_builds_rows_and_aggregate_counts() -> None:
    uploads = [
        make_upload("a.ab1"),
        make_upload("b.ab1"),
        make_upload("broken.ab1"),
    ]
    trace_data = TraceData(channels={"DATA9": [1, 2, 3]})

    a_result = trim_sequence_record(
        make_record(
            name="trace_a",
            sequence="ACGTACGT",
            qualities=[5, 8, 30, 35, 40, 25, 7, 3],
            trace_data=trace_data,
        ),
        TrimConfig(
            left_trim=1,
            right_trim=1,
            quality_trim_enabled=True,
            error_probability_cutoff=0.01,
        ),
    )
    b_result = trim_sequence_record(
        make_record(
            name="trace_b",
            sequence="ACGT",
            qualities=None,
        ),
        TrimConfig(left_trim=1, right_trim=1, min_length=5),
    )

    summary = build_batch_summary(
        uploads=uploads,
        trim_results={
            "a.ab1": a_result,
            "b.ab1": b_result,
        },
        parse_errors={"broken.ab1": "Failed to parse ABI file: broken.ab1"},
    )

    assert summary.total_uploaded_files == 3
    assert summary.parsed_files == 2
    assert summary.failed_files == 1
    assert summary.trimmed_records == 2
    assert summary.records_passing_min_length == 1
    assert summary.records_failing_min_length == 1
    assert summary.fastq_exportable_records == 1
    assert summary.parse_errors == (
        ("broken.ab1", "Failed to parse ABI file: broken.ab1"),
    )

    a_row, b_row = summary.records

    assert a_row.source_filename == "a.ab1"
    assert a_row.record_id == "trace_a_id"
    assert a_row.display_name == "trace_a"
    assert a_row.original_length == 8
    assert a_row.trimmed_length == 2
    assert a_row.orientation == "forward"
    assert a_row.total_bases_removed == 6
    assert a_row.fixed_bases_removed_left == 1
    assert a_row.fixed_bases_removed_right == 1
    assert a_row.quality_bases_removed_left == 2
    assert a_row.quality_bases_removed_right == 2
    assert a_row.passed_min_length is True
    assert a_row.has_qualities is True
    assert a_row.has_trace_data is True
    assert a_row.fastq_exportable is True
    assert a_row.to_row()["orientation"] == "forward"

    assert b_row.source_filename == "b.ab1"
    assert b_row.trimmed_length == 2
    assert b_row.orientation == "forward"
    assert b_row.total_bases_removed == 2
    assert b_row.fixed_bases_removed_left == 1
    assert b_row.fixed_bases_removed_right == 1
    assert b_row.quality_bases_removed_left == 0
    assert b_row.quality_bases_removed_right == 0
    assert b_row.passed_min_length is False
    assert b_row.has_qualities is False
    assert b_row.has_trace_data is False
    assert b_row.fastq_exportable is False


def test_build_batch_summary_handles_all_failures() -> None:
    uploads = [make_upload("broken.ab1")]

    summary = build_batch_summary(
        uploads=uploads,
        trim_results={},
        parse_errors={"broken.ab1": "Failed to parse ABI file: broken.ab1"},
    )

    assert summary.total_uploaded_files == 1
    assert summary.parsed_files == 0
    assert summary.failed_files == 1
    assert summary.trimmed_records == 0
    assert summary.records_passing_min_length == 0
    assert summary.records_failing_min_length == 0
    assert summary.fastq_exportable_records == 0
    assert summary.records == ()
    assert summary.table_rows() == []


def test_build_batch_export_policy_partitions_records_and_reasons() -> None:
    a_result = trim_sequence_record(
        make_record(
            name="trace_a",
            sequence="ACGT",
            qualities=[40, 41, 42, 43],
        ),
        TrimConfig(),
    )
    b_result = trim_sequence_record(
        make_record(
            name="trace_b",
            sequence="ACGT",
            qualities=None,
        ),
        TrimConfig(min_length=5),
    )

    policy = build_batch_export_policy(
        trim_results={
            "a.ab1": a_result,
            "b.ab1": b_result,
        }
    )

    assert tuple(record.source_filename for record in policy.fasta_records) == (
        "a.ab1",
        "b.ab1",
    )
    assert tuple(record.source_filename for record in policy.fastq_records) == (
        "a.ab1",
    )
    assert tuple(
        record.source_filename for record in policy.passing_min_length_records
    ) == ("a.ab1",)
    assert tuple(
        record.source_filename for record in policy.failing_min_length_records
    ) == ("b.ab1",)

    assert policy.eligible_records(export_format="fastq") == (a_result.record,)
    assert policy.eligible_records(
        export_format="fasta",
        require_min_length=True,
    ) == (a_result.record,)

    assert policy.ineligible_reasons_by_filename(export_format="fastq") == (
        ("b.ab1", ("missing per-base qualities",)),
    )
    assert policy.ineligible_reasons_by_filename(
        export_format="fasta",
        require_min_length=True,
    ) == (("b.ab1", ("below minimum length",)),)
    assert policy.ineligible_reasons_by_filename(
        export_format="fastq",
        require_min_length=True,
    ) == (("b.ab1", ("missing per-base qualities", "below minimum length")),)


def test_build_batch_export_policy_marks_misaligned_fastq_records_ineligible() -> None:
    bad_record = make_record(
        name="trace_bad",
        sequence="ACGT",
        qualities=[40, 41, 42],
    )
    bad_result = TrimResult(
        record=bad_record,
        original_length=4,
        trimmed_length=4,
        bases_removed_left=0,
        bases_removed_right=0,
        passed_min_length=True,
    )

    policy = build_batch_export_policy(trim_results={"bad.ab1": bad_result})

    assert tuple(record.source_filename for record in policy.fastq_records) == ()
    assert policy.eligible_records(export_format="fasta") == (bad_record,)
    assert policy.ineligible_reasons_by_filename(export_format="fastq") == (
        ("bad.ab1", ("sequence/quality length mismatch",)),
    )
