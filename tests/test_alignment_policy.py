from __future__ import annotations

from abi_sauce.alignment_policy import (
    alignment_overlap_metrics,
    alignment_strands_for_policy,
    build_semiglobal_aligner,
    oriented_qualities,
    select_best_oriented_alignment,
)


def test_alignment_strands_for_policy() -> None:
    assert alignment_strands_for_policy("auto") == ("forward", "reverse_complement")
    assert alignment_strands_for_policy("forward") == ("forward",)
    assert alignment_strands_for_policy("reverse_complement") == ("reverse_complement",)


def test_oriented_qualities_reverses_reverse_complement() -> None:
    assert oriented_qualities([10, 20, 30], strand="forward") == [10, 20, 30]
    assert oriented_qualities([10, 20, 30], strand="reverse_complement") == [30, 20, 10]
    assert oriented_qualities(None, strand="forward") is None


def test_select_best_oriented_alignment_prefers_reverse_complement_when_needed() -> (
    None
):
    aligner = build_semiglobal_aligner()
    result = select_best_oriented_alignment(
        target_sequence="CCCCAAAAC",
        display_query_sequence="GTTTT",
        display_query_qualities=[40] * 5,
        aligner=aligner,
    )

    assert result is not None
    assert result.strand == "reverse_complement"
    assert result.sequence == "AAAAC"
    assert result.qualities == [40] * 5


def test_alignment_overlap_metrics_reports_overlap_and_identity() -> None:
    aligner = build_semiglobal_aligner()
    alignment = aligner.align("AACCGG", "AACCTG")[0]

    overlap_length, percent_identity = alignment_overlap_metrics(
        alignment,
        target_sequence="AACCGG",
        query_sequence="AACCTG",
    )

    assert overlap_length == 6
    assert percent_identity == (5 / 6) * 100.0
