from __future__ import annotations

from abi_sauce.alignment_state import AlignmentDefinition
from abi_sauce.assembly_types import AssemblyConfig, AssemblyResult
from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData
from abi_sauce.services.alignment_compute import (
    compute_saved_alignment,
    compute_saved_alignments,
)
from abi_sauce.services.batch_parse import ParsedBatch, build_batch_signature
from abi_sauce.services.batch_trim import apply_trim_configs
from abi_sauce.trimming import TrimConfig


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
        description="synthetic alignment record",
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


def make_prepared_batch():
    uploads = (
        SequenceUpload(filename="left.ab1", content=b"left"),
        SequenceUpload(filename="right.ab1", content=b"right"),
        SequenceUpload(filename="trace.ab1", content=b"trace"),
        SequenceUpload(filename="trace_2.ab1", content=b"trace_2"),
        SequenceUpload(filename="trace_3.ab1", content=b"trace_3"),
    )
    parsed_batch = ParsedBatch(
        uploads=uploads,
        parsed_records={
            "left.ab1": make_record(
                name="left",
                sequence="CCCCAAAAC",
                qualities=[40] * 9,
            ),
            "right.ab1": make_record(
                name="right",
                sequence="GTTTT",
                qualities=[40] * 5,
            ),
            "trace.ab1": make_record(
                name="trace",
                sequence="AACCGGTT",
                qualities=[10, 20, 30, 40, 50, 40, 30, 20],
            ),
            "trace_2.ab1": make_record(
                name="trace_2",
                sequence="AACCGGTT",
                qualities=[40] * 8,
            ),
            "trace_3.ab1": make_record(
                name="trace_3",
                sequence="AACCGGTA",
                qualities=[35] * 8,
            ),
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )
    return apply_trim_configs(
        parsed_batch,
        trim_configs_by_name={"trace.ab1": TrimConfig(left_trim=2, right_trim=2)},
    )


def test_compute_saved_alignment_returns_pairwise_assembly_result() -> None:
    prepared_batch = make_prepared_batch()
    definition = AlignmentDefinition(
        alignment_id="alignment-pairwise",
        name="Amplicon A",
        source_filenames=("left.ab1", "right.ab1"),
        engine_kind="pairwise",
        assembly_config=AssemblyConfig(
            min_overlap_length=4,
            min_percent_identity=90.0,
        ),
    )

    computed_alignment = compute_saved_alignment(prepared_batch, definition)

    assert computed_alignment.status == "ok"
    assert computed_alignment.assembly is not None
    assert computed_alignment.assembly.consensus_record is not None
    assert computed_alignment.assembly.consensus_record.sequence == "CCCCAAAAC"
    assert isinstance(computed_alignment.assembly.result, AssemblyResult)


def test_compute_saved_alignment_returns_reference_alignment_result() -> None:
    prepared_batch = make_prepared_batch()
    definition = AlignmentDefinition(
        alignment_id="alignment-reference",
        name="Trace vs ref",
        source_filenames=("trace.ab1",),
        engine_kind="reference_single",
        reference_name="ref",
        reference_text=">ref\nCCGA\n",
        strand_policy="forward",
    )

    computed_alignment = compute_saved_alignment(prepared_batch, definition)

    assert computed_alignment.status == "ok"
    assert computed_alignment.reference_alignment is not None
    assert (
        computed_alignment.reference_alignment.alignment_result.reference_name == "ref"
    )
    assert computed_alignment.reference_alignment.alignment_result.mismatch_count == 1
    assert len(computed_alignment.reference_alignment.event_rows) == 1


def test_compute_saved_alignment_rejects_missing_reference_text() -> None:
    prepared_batch = make_prepared_batch()
    definition = AlignmentDefinition(
        alignment_id="alignment-reference",
        name="Trace vs ref",
        source_filenames=("trace.ab1",),
        engine_kind="reference_single",
        reference_text=None,
    )

    computed_alignment = compute_saved_alignment(prepared_batch, definition)

    assert computed_alignment.status == "invalid_definition"
    assert computed_alignment.reference_alignment is None
    assert computed_alignment.status_reason is not None
    assert "requires reference text" in computed_alignment.status_reason


def test_compute_saved_alignment_returns_reference_multi_alignment_result() -> None:
    prepared_batch = make_prepared_batch()
    definition = AlignmentDefinition(
        alignment_id="alignment-reference-multi",
        name="Many vs ref",
        source_filenames=("trace_2.ab1", "trace_3.ab1"),
        engine_kind="reference_multi",
        reference_name="ref",
        reference_text=">ref\nAACCGGTT\n",
        strand_policy="forward",
        assembly_config=AssemblyConfig(
            min_overlap_length=6,
            min_percent_identity=80.0,
            quality_margin=3,
        ),
    )

    computed_alignment = compute_saved_alignment(prepared_batch, definition)

    assert computed_alignment.status == "ok"
    assert computed_alignment.reference_multi_alignment is not None
    assert computed_alignment.reference_multi_alignment.result.reference_name == "ref"
    assert (
        computed_alignment.reference_multi_alignment.result.included_member_count == 2
    )
    assert (
        computed_alignment.reference_multi_alignment.result.consensus_sequence
        == "AACCGGTT"
    )


def test_compute_saved_alignments_preserves_definition_order() -> None:
    prepared_batch = make_prepared_batch()
    definitions = (
        AlignmentDefinition(
            alignment_id="alignment-pairwise",
            name="Amplicon A",
            source_filenames=("left.ab1", "right.ab1"),
            engine_kind="pairwise",
            assembly_config=AssemblyConfig(
                min_overlap_length=4,
                min_percent_identity=90.0,
            ),
        ),
        AlignmentDefinition(
            alignment_id="alignment-reference",
            name="Trace vs ref",
            source_filenames=("trace.ab1",),
            engine_kind="reference_single",
            reference_name="ref",
            reference_text=">ref\nCCGA\n",
            strand_policy="forward",
        ),
    )

    computed_alignments = compute_saved_alignments(prepared_batch, definitions)

    assert tuple(computed_alignments) == (
        "alignment-pairwise",
        "alignment-reference",
    )
