# tests/test_reference_alignment.py
from __future__ import annotations

from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment import (
    align_trimmed_read_to_reference,
    alignment_events_to_rows,
    normalize_reference,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    trace_length = 100
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=TraceData(
            channels={
                "DATA9": [1] * trace_length,
                "DATA10": [1] * trace_length,
                "DATA11": [1] * trace_length,
                "DATA12": [1] * trace_length,
            },
            base_positions=base_positions
            or list(range(5, 5 + (10 * len(sequence)), 10)),
            channel_order="GATC",
        ),
    )


def test_normalize_reference_parses_fasta_header_and_sequence() -> None:
    reference_name, reference_sequence = normalize_reference(
        ">amplicon_001\nACGT AC\nGT\n"
    )

    assert reference_name == "amplicon_001"
    assert reference_sequence == "ACGTACGT"


def test_align_trimmed_read_to_reference_returns_perfect_forward_match() -> None:
    raw_record = make_record(
        name="trace_forward",
        sequence="AACCGGTT",
        qualities=[10, 20, 30, 40, 50, 40, 30, 20],
    )
    trim_result = trim_sequence_record(
        raw_record,
        TrimConfig(left_trim=2, right_trim=2),
    )

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=">ref\nCCGG\n",
        strand_policy="auto",
    )

    assert alignment_result.strand == "forward"
    assert alignment_result.reference_start == 1
    assert alignment_result.reference_end == 4
    assert alignment_result.query_start == 1
    assert alignment_result.query_end == 4
    assert alignment_result.percent_identity == 100.0
    assert alignment_result.mismatch_count == 0
    assert alignment_result.insertion_count == 0
    assert alignment_result.deletion_count == 0
    assert alignment_events_to_rows(alignment_result) == []


def test_align_trimmed_read_to_reference_auto_picks_reverse_complement() -> None:
    raw_record = make_record(
        name="trace_reverse",
        sequence="GTTTT",
        qualities=[10, 11, 12, 13, 14],
    )
    trim_result = trim_sequence_record(raw_record, TrimConfig())

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text="AAAAC",
        strand_policy="auto",
    )

    assert alignment_result.strand == "reverse-complement"
    assert alignment_result.percent_identity == 100.0
    assert alignment_result.aligned_reference == "AAAAC"
    assert alignment_result.aligned_query == "AAAAC"


def test_align_trimmed_read_to_reference_forward_uses_display_orientation() -> None:
    raw_record = make_record(
        name="trace_display_forward",
        sequence="GTTTT",
        qualities=[10, 11, 12, 13, 14],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(raw_record, TrimConfig())

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text="AAAAC",
        strand_policy="forward",
    )

    assert alignment_result.strand == "forward"
    assert alignment_result.percent_identity == 100.0
    assert alignment_result.aligned_reference == "AAAAC"
    assert alignment_result.aligned_query == "AAAAC"


def test_alignment_events_include_qscore_and_trace_position_for_mismatch() -> None:
    raw_record = make_record(
        name="trace_mismatch",
        sequence="AACCGGTT",
        qualities=[10, 20, 30, 40, 50, 40, 30, 20],
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    trim_result = trim_sequence_record(
        raw_record,
        TrimConfig(left_trim=2, right_trim=2),
    )

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text="CCGA",
        strand_policy="forward",
    )
    event_rows = alignment_events_to_rows(alignment_result)

    assert len(event_rows) == 1
    assert event_rows[0]["type"] == "mismatch"
    assert event_rows[0]["ref_pos"] == 4
    assert event_rows[0]["query_pos"] == 4
    assert event_rows[0]["ref_base"] == "A"
    assert event_rows[0]["query_base"] == "G"
    assert event_rows[0]["qscore"] == 40
    assert event_rows[0]["trace_x"] == 55


def test_alignment_events_use_display_trace_coordinates_for_reverse_complement_samples() -> (
    None
):
    raw_record = make_record(
        name="trace_display_trace_x",
        sequence="GTTTT",
        qualities=[10, 20, 30, 40, 50],
        base_positions=[5, 15, 25, 35, 45],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(raw_record, TrimConfig())

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text="AAAAT",
        strand_policy="forward",
    )
    event_rows = alignment_events_to_rows(alignment_result)

    assert len(event_rows) == 1
    assert event_rows[0]["type"] == "mismatch"
    assert event_rows[0]["query_pos"] == 5
    assert event_rows[0]["trace_x"] == 94
