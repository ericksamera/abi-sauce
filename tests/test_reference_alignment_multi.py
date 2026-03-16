from __future__ import annotations

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment_multi import align_trimmed_reads_to_reference
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    resolved_base_positions = base_positions or list(
        range(5, 5 + (10 * len(sequence)), 10)
    )
    trace_length = max(resolved_base_positions[-1] + 10, 200)
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic shared-reference alignment record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=TraceData(
            channels={
                "DATA9": list(range(trace_length)),
                "DATA10": [value * 2 for value in range(trace_length)],
                "DATA11": [value * 3 for value in range(trace_length)],
                "DATA12": [value * 4 for value in range(trace_length)],
            },
            base_positions=resolved_base_positions,
            channel_order="GATC",
        ),
    )


def make_inputs():
    raw_records_by_source_filename = {
        "read_1.ab1": make_record(
            name="read_1",
            sequence="AACCGGTT",
            qualities=[40] * 8,
        ),
        "read_2.ab1": make_record(
            name="read_2",
            sequence="AACCGGTT",
            qualities=[38] * 8,
        ),
        "read_3.ab1": make_record(
            name="read_3",
            sequence="AACCGGTA",
            qualities=[35] * 8,
        ),
        "read_ins.ab1": make_record(
            name="read_ins",
            sequence="AACCGGGTT",
            qualities=[30] * 9,
        ),
        "read_bad.ab1": make_record(
            name="read_bad",
            sequence="TTTTTTTT",
            qualities=[20] * 8,
        ),
    }
    trim_results_by_source_filename = {
        source_filename: trim_sequence_record(record, TrimConfig())
        for source_filename, record in raw_records_by_source_filename.items()
    }
    return raw_records_by_source_filename, trim_results_by_source_filename


def test_align_trimmed_reads_to_reference_builds_shared_reference_consensus() -> None:
    raw_records_by_source_filename, trim_results_by_source_filename = make_inputs()

    result = align_trimmed_reads_to_reference(
        source_filenames=("read_1.ab1", "read_2.ab1", "read_3.ab1", "read_bad.ab1"),
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
        reference_text=">ref\nAACCGGTT\n",
        strand_policy="forward",
        config=AssemblyConfig(
            min_overlap_length=6,
            min_percent_identity=80.0,
            quality_margin=3,
        ),
    )

    assert result.reference_name == "ref"
    assert result.accepted is True
    assert result.included_member_count == 3
    assert result.excluded_member_count == 1
    assert result.gapped_reference == "AACCGGTT"
    assert result.consensus_sequence == "AACCGGTT"
    assert result.columns[-1].ref_base == "T"
    assert result.columns[-1].support_counts == (("T", 2), ("A", 1))
    assert result.columns[-1].consensus_base == "T"
    assert result.columns[-1].resolution == "majority_resolved"

    excluded_member = next(member for member in result.members if not member.included)
    assert excluded_member.source_filename == "read_bad.ab1"
    assert excluded_member.inclusion_reason is not None
    assert excluded_member.inclusion_reason.startswith("overlap length below threshold")


def test_align_trimmed_reads_to_reference_resolves_insertion_buckets() -> None:
    raw_records_by_source_filename, trim_results_by_source_filename = make_inputs()

    result = align_trimmed_reads_to_reference(
        source_filenames=("read_1.ab1", "read_2.ab1", "read_ins.ab1"),
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
        reference_text=">ref\nAACCGGTT\n",
        strand_policy="forward",
        config=AssemblyConfig(
            min_overlap_length=6,
            min_percent_identity=80.0,
            quality_margin=3,
        ),
    )

    insertion_column = next(
        column for column in result.columns if column.anchor_kind == "insertion"
    )

    assert result.accepted is True
    assert result.gapped_reference.count("-") == 1
    assert insertion_column.ref_base == "-"
    assert insertion_column.consensus_base != "-"
    assert insertion_column.resolution == "single_read"
    assert insertion_column.non_gap_member_count == 1
    assert insertion_column.gap_member_count == 2
