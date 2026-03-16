from __future__ import annotations

from abi_sauce.models import SequenceRecord
from abi_sauce.oriented_reads import (
    nonempty_prepared_reads,
    prepare_trimmed_read,
    prepare_trimmed_reads,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        orientation="forward",
        trace_data=None,
        annotations={"source_filename": f"{name}.ab1"},
    )


def test_prepare_trimmed_read_materializes_display_orientation() -> None:
    raw_record = make_record(
        name="sample",
        sequence="GTTTT",
        qualities=[10, 20, 30, 40, 50],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(raw_record, TrimConfig())

    prepared_read = prepare_trimmed_read(
        raw_record=raw_record,
        trim_result=trim_result,
    )

    assert prepared_read.source_filename == "sample.ab1"
    assert prepared_read.display_name == "sample"
    assert prepared_read.display_sequence == "AAAAC"
    assert prepared_read.display_qualities == [50, 40, 30, 20, 10]
    assert prepared_read.trimmed_length == 5
    assert prepared_read.has_trace_data is False
    assert prepared_read.has_qualities is True


def test_prepare_trimmed_reads_preserves_selection_order() -> None:
    left_record = make_record(name="left", sequence="AACCGG")
    right_record = make_record(name="right", sequence="TTGGCC")
    prepared_reads = prepare_trimmed_reads(
        source_filenames=("right.ab1", "left.ab1"),
        raw_records_by_source_filename={
            "left.ab1": left_record,
            "right.ab1": right_record,
        },
        trim_results_by_source_filename={
            "left.ab1": trim_sequence_record(left_record, TrimConfig()),
            "right.ab1": trim_sequence_record(right_record, TrimConfig()),
        },
    )

    assert tuple(prepared_read.source_filename for prepared_read in prepared_reads) == (
        "right.ab1",
        "left.ab1",
    )


def test_nonempty_prepared_reads_filters_empty_display_sequences() -> None:
    kept_record = make_record(name="kept", sequence="AACCGG")
    empty_record = make_record(name="empty", sequence="AACCGG")
    empty_trim_result = trim_sequence_record(empty_record, TrimConfig(left_trim=6))

    prepared_reads = (
        prepare_trimmed_read(
            source_filename="kept.ab1",
            raw_record=kept_record,
            trim_result=trim_sequence_record(kept_record, TrimConfig()),
        ),
        prepare_trimmed_read(
            source_filename="empty.ab1",
            raw_record=empty_record,
            trim_result=empty_trim_result,
        ),
    )

    assert tuple(
        prepared_read.source_filename
        for prepared_read in nonempty_prepared_reads(prepared_reads)
    ) == ("kept.ab1",)
