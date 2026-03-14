from __future__ import annotations

from abi_sauce.assembly import (
    AssemblyConfig,
    assemble_trimmed_pair,
    assembly_conflicts_to_rows,
    consensus_record_from_result,
)
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    trace_length = 200
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic assembly record",
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


def test_assemble_trimmed_pair_auto_picks_reverse_complement_for_right_read() -> None:
    left_raw_record = make_record(
        name="left",
        sequence="CCCCAAAAC",
        qualities=[40] * 9,
    )
    right_raw_record = make_record(
        name="right",
        sequence="GTTTT",
        qualities=[40] * 5,
    )

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig())
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
    )

    assert result.accepted is True
    assert result.chosen_right_orientation == "reverse-complement"
    assert result.overlap_length == 5
    assert result.percent_identity == 100.0
    assert result.consensus_sequence == "CCCCAAAAC"


def test_assemble_trimmed_pair_rejects_insufficient_overlap() -> None:
    left_raw_record = make_record(name="left", sequence="AAAACCCC", qualities=[40] * 8)
    right_raw_record = make_record(
        name="right", sequence="CCCCGGGG", qualities=[40] * 8
    )

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig())
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=5, min_percent_identity=90.0),
    )

    assert result.accepted is False
    assert result.rejection_reason is not None
    assert "overlap length below threshold" in result.rejection_reason


def test_assemble_trimmed_pair_resolves_mismatch_by_quality() -> None:
    left_raw_record = make_record(
        name="left",
        sequence="AAAACCCCTTTT",
        qualities=[40] * 9 + [35] + [40] * 2,
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115],
    )
    right_raw_record = make_record(
        name="right",
        sequence="CCCCTGTT",
        qualities=[40, 40, 40, 40, 40, 10, 40, 40],
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig(left_trim=4))
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=8, min_percent_identity=70.0),
    )

    assert result.accepted is True
    assert result.consensus_sequence == "CCCCTTTT"
    assert result.conflict_count == 1
    conflict = result.conflicts[0]
    assert conflict.resolution == "quality_resolved"
    assert conflict.consensus_base == "T"
    assert conflict.left_quality == 35
    assert conflict.right_quality == 10
    assert conflict.left_trace_x == 95
    assert conflict.right_trace_x == 55


def test_assemble_trimmed_pair_emits_N_for_unresolved_mismatch_without_qualities() -> (
    None
):
    left_raw_record = make_record(name="left", sequence="AAAACCCCTTTT", qualities=None)
    right_raw_record = make_record(name="right", sequence="CCCCTGTT", qualities=None)

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig(left_trim=4))
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=8, min_percent_identity=70.0),
    )

    assert result.accepted is True
    assert result.consensus_sequence == "CCCCTNTT"
    assert result.conflict_count == 1
    conflict = result.conflicts[0]
    assert conflict.resolution == "ambiguous"
    assert conflict.consensus_base == "N"


def test_consensus_record_projection_and_conflict_rows() -> None:
    left_raw_record = make_record(name="left", sequence="AAAACCCCTTTT", qualities=None)
    right_raw_record = make_record(name="right", sequence="CCCCTGTT", qualities=None)

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig(left_trim=4))
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=8, min_percent_identity=70.0),
    )

    consensus_record = consensus_record_from_result(result)
    rows = assembly_conflicts_to_rows(result)

    assert consensus_record.sequence == "CCCCTNTT"
    assert consensus_record.source_format == "assembly"
    assert consensus_record.annotations["assembly_conflict_count"] == 1
    assert consensus_record.annotations["assembly_accepted"] is True
    assert rows[0]["consensus_base"] == "N"
    assert rows[0]["resolution"] == "ambiguous"


def test_consensus_record_projection_accepts_explicit_name() -> None:
    left_raw_record = make_record(name="left", sequence="AAAACCCCTTTT", qualities=None)
    right_raw_record = make_record(name="right", sequence="CCCCTGTT", qualities=None)

    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig(left_trim=4))
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=8, min_percent_identity=70.0),
    )

    consensus_record = consensus_record_from_result(result, name="Amplicon 1")

    assert consensus_record.name == "Amplicon 1"
    assert consensus_record.record_id == "Amplicon 1"
