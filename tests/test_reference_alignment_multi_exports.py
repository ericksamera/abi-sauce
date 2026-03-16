from __future__ import annotations

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment_multi import align_trimmed_reads_to_reference
from abi_sauce.reference_alignment_multi_exports import (
    consensus_record_from_reference_multi_result,
    format_reference_multi_alignment_fasta,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(name: str, sequence: str) -> SequenceRecord:
    trace_length = 200
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic shared-reference export record",
        sequence=sequence,
        source_format="abi",
        qualities=[40] * len(sequence),
        trace_data=TraceData(
            channels={
                "DATA9": [1] * trace_length,
                "DATA10": [1] * trace_length,
                "DATA11": [1] * trace_length,
                "DATA12": [1] * trace_length,
            },
            base_positions=list(range(5, 5 + (10 * len(sequence)), 10)),
            channel_order="GATC",
        ),
    )


def make_result():
    raw_records_by_source_filename = {
        "read_1.ab1": make_record("read_1", "AACCGGTT"),
        "read_2.ab1": make_record("read_2", "AACCGGTT"),
        "read_3.ab1": make_record("read_3", "AACCGGTA"),
    }
    trim_results_by_source_filename = {
        source_filename: trim_sequence_record(record, TrimConfig())
        for source_filename, record in raw_records_by_source_filename.items()
    }
    return align_trimmed_reads_to_reference(
        source_filenames=("read_1.ab1", "read_2.ab1", "read_3.ab1"),
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
        reference_text=">ref\nAACCGGTT\n",
        strand_policy="forward",
        config=AssemblyConfig(min_overlap_length=6, min_percent_identity=80.0),
    )


def test_format_reference_multi_alignment_fasta_renders_reference_members_and_consensus() -> (
    None
):
    result = make_result()

    fasta = format_reference_multi_alignment_fasta(result, consensus_name="consensus")

    assert fasta.startswith(">ref\nAACCGGTT\n>read_1__forward\nAACCGGTT\n")
    assert ">consensus__gapped_consensus\nAACCGGTT\n" in fasta


def test_consensus_record_from_reference_multi_result_projects_consensus_sequence() -> (
    None
):
    result = make_result()

    consensus_record = consensus_record_from_reference_multi_result(
        result,
        name="shared_consensus",
    )

    assert consensus_record.name == "shared_consensus"
    assert consensus_record.sequence == "AACCGGTT"
    assert consensus_record.annotations["alignment_engine_kind"] == "reference_multi"
