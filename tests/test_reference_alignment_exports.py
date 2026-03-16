from __future__ import annotations

from abi_sauce.reference_alignment_exports import format_reference_alignment_fasta
from abi_sauce.reference_alignment_types import AlignmentEvent, AlignmentResult


def test_format_reference_alignment_fasta_renders_gapped_entries() -> None:
    result = AlignmentResult(
        sample_name="trace",
        reference_name="amplicon_001",
        strand="forward",
        score=10.0,
        reference_start=1,
        reference_end=4,
        query_start=1,
        query_end=4,
        percent_identity=75.0,
        mismatch_count=1,
        insertion_count=0,
        deletion_count=0,
        aligned_reference="A-CGT",
        match_line="| .||",
        aligned_query="ATCGT",
        events=(
            AlignmentEvent(
                ref_pos=2,
                query_pos=2,
                event_type="mismatch",
                ref_base="C",
                query_base="T",
                qscore=30,
                flank_q_left=20,
                flank_q_right=40,
                trace_x=25,
                context_ref="ACGT",
                context_query="ATGT",
            ),
        ),
    )

    fasta = format_reference_alignment_fasta(result)

    assert fasta == (">amplicon_001\n" "A-CGT\n" ">trace__forward\n" "ATCGT\n")
