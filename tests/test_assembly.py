from __future__ import annotations

from abi_sauce.assembly import (
    AssemblyConfig,
    assemble_trimmed_multi,
    assemble_trimmed_pair,
    assembly_conflicts_to_rows,
    consensus_record_from_multi_result,
    consensus_record_from_result,
    format_assembly_alignment_fasta,
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
    assert result.chosen_right_orientation == "reverse_complement"
    assert result.overlap_length == 5
    assert result.percent_identity == 100.0
    assert result.consensus_sequence == "CCCCAAAAC"
    assert len(result.columns) == len(result.aligned_left)
    assert tuple(column.column_index for column in result.columns) == tuple(
        range(1, len(result.aligned_left) + 1)
    )
    assert tuple(column.consensus_base for column in result.columns) == tuple(
        result.gapped_consensus
    )


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
    assert len(result.columns) == 8
    mismatch_column = result.columns[5]
    assert mismatch_column.column_index == 6
    assert mismatch_column.left_query_index == 5
    assert mismatch_column.right_query_index == 5
    assert mismatch_column.left_query_pos == 6
    assert mismatch_column.right_query_pos == 6
    assert mismatch_column.left_base == "T"
    assert mismatch_column.right_base == "G"
    assert mismatch_column.consensus_base == "T"
    assert mismatch_column.resolution == "quality_resolved"
    assert mismatch_column.is_overlap is True
    assert mismatch_column.is_match is False
    assert mismatch_column.left_trace_x == 95
    assert mismatch_column.right_trace_x == 55
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


def test_assemble_trimmed_multi_builds_seeded_consensus_and_alignment_grid() -> None:
    seed_raw_record = make_record(
        name="seed",
        sequence="AACCGGTTA",
        qualities=[40] * 9,
    )
    shifted_raw_record = make_record(
        name="shifted",
        sequence="ACCGGATTA",
        qualities=[35] * 9,
    )
    reverse_raw_record = make_record(
        name="reverse",
        sequence="TAACCGGTT",
        qualities=[30] * 9,
    )

    trim_results = {
        "seed.ab1": trim_sequence_record(seed_raw_record, TrimConfig()),
        "shifted.ab1": trim_sequence_record(shifted_raw_record, TrimConfig()),
        "reverse.ab1": trim_sequence_record(reverse_raw_record, TrimConfig()),
    }
    raw_records = {
        "seed.ab1": seed_raw_record,
        "shifted.ab1": shifted_raw_record,
        "reverse.ab1": reverse_raw_record,
    }

    result = assemble_trimmed_multi(
        source_filenames=("seed.ab1", "shifted.ab1", "reverse.ab1"),
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=70.0),
    )

    assert result.accepted is True
    assert result.seed_member_index == 0
    assert result.included_member_indices == (0, 1, 2)
    assert result.included_member_count == 3
    assert result.excluded_member_count == 0
    assert result.gapped_consensus == "AACCGGATTA"
    assert result.aligned_member_sequences == (
        "AACCGG-TTA",
        "-ACCGGATTA",
        "AACCGG-TTA",
    )
    assert result.members[0].is_seed is True
    assert result.members[2].chosen_orientation == "reverse_complement"
    assert len(result.columns) == 10
    assert result.columns[6].consensus_base == "A"
    assert result.columns[6].resolution == "single_read"

    consensus_record = consensus_record_from_multi_result(result, name="Amplicon M")
    assert consensus_record.name == "Amplicon M"
    assert consensus_record.sequence == "AACCGGATTA"
    assert consensus_record.annotations["assembly_engine_kind"] == "multi"
    assert consensus_record.annotations["assembly_included_member_count"] == 3


def test_assemble_trimmed_multi_excludes_low_identity_members() -> None:
    seed_raw_record = make_record(
        name="seed",
        sequence="AACCGGTTA",
        qualities=[40] * 9,
    )
    reverse_raw_record = make_record(
        name="reverse",
        sequence="TAACCGGTT",
        qualities=[30] * 9,
    )
    bad_raw_record = make_record(
        name="bad",
        sequence="GGGGGGGGG",
        qualities=[20] * 9,
    )

    trim_results = {
        "seed.ab1": trim_sequence_record(seed_raw_record, TrimConfig()),
        "reverse.ab1": trim_sequence_record(reverse_raw_record, TrimConfig()),
        "bad.ab1": trim_sequence_record(bad_raw_record, TrimConfig()),
    }
    raw_records = {
        "seed.ab1": seed_raw_record,
        "reverse.ab1": reverse_raw_record,
        "bad.ab1": bad_raw_record,
    }

    result = assemble_trimmed_multi(
        source_filenames=("seed.ab1", "reverse.ab1", "bad.ab1"),
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=70.0),
    )

    assert result.accepted is True
    assert result.included_member_count == 2
    assert result.excluded_member_count == 1
    assert result.consensus_sequence == "AACCGGTTA"
    assert result.members[2].included is False
    assert (
        result.members[2].inclusion_reason
        == "no overlapping aligned columns were found"
    )


def test_format_assembly_alignment_fasta_includes_pairwise_rows_and_gapped_consensus() -> (
    None
):
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

    alignment_fasta = format_assembly_alignment_fasta(
        result,
        consensus_name="Amplicon A",
    )
    lines = alignment_fasta.strip().splitlines()

    assert lines == [
        f">{result.left_display_name}",
        result.aligned_left,
        f">{result.right_display_name}",
        result.aligned_right,
        ">Amplicon A__gapped_consensus",
        result.gapped_consensus,
    ]


def test_format_assembly_alignment_fasta_includes_multi_rows_and_gapped_consensus() -> (
    None
):
    seed_raw_record = make_record(
        name="seed",
        sequence="AACCGGTTA",
        qualities=[40] * 9,
    )
    shifted_raw_record = make_record(
        name="shifted",
        sequence="ACCGGATTA",
        qualities=[35] * 9,
    )
    reverse_raw_record = make_record(
        name="reverse",
        sequence="TAACCGGTT",
        qualities=[30] * 9,
    )

    trim_results = {
        "seed.ab1": trim_sequence_record(seed_raw_record, TrimConfig()),
        "shifted.ab1": trim_sequence_record(shifted_raw_record, TrimConfig()),
        "reverse.ab1": trim_sequence_record(reverse_raw_record, TrimConfig()),
    }
    raw_records = {
        "seed.ab1": seed_raw_record,
        "shifted.ab1": shifted_raw_record,
        "reverse.ab1": reverse_raw_record,
    }

    result = assemble_trimmed_multi(
        source_filenames=("seed.ab1", "shifted.ab1", "reverse.ab1"),
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=70.0),
    )

    alignment_fasta = format_assembly_alignment_fasta(
        result,
        consensus_name="Amplicon M",
    )
    lines = alignment_fasta.strip().splitlines()

    assert lines == [
        ">seed",
        result.aligned_member_sequences[0],
        ">shifted",
        result.aligned_member_sequences[1],
        ">reverse",
        result.aligned_member_sequences[2],
        ">Amplicon M__gapped_consensus",
        result.gapped_consensus,
    ]
