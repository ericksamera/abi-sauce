from __future__ import annotations

from abi_sauce.models import SequenceOrientation, SequenceRecord
from abi_sauce.orientation import (
    complement_base,
    complement_sequence,
    materialize_oriented_record,
    orient_left_right_values,
    orient_qualities,
    orient_sequence,
    orient_trim_config_for_display,
    raw_trim_config_from_display,
    reverse_complement_sequence,
)
from abi_sauce.trimming import TrimConfig


def make_record(
    *,
    sequence: str = "ACGT",
    qualities: list[int] | None = None,
    orientation: SequenceOrientation = "forward",
) -> SequenceRecord:
    return SequenceRecord(
        record_id="trace_001",
        name="trace_001",
        description="test record",
        sequence=sequence,
        source_format="abi",
        orientation=orientation,
        qualities=qualities,
    )


def test_complement_sequence_supports_iupac_bases_without_reversing() -> None:
    assert complement_sequence("ACGTRYKMSWBDHVN") == "TGCAYRMKSWVHDBN"


def test_complement_base_maps_single_base_labels() -> None:
    assert complement_base("A") == "T"


def test_reverse_complement_sequence_supports_iupac_bases() -> None:
    assert reverse_complement_sequence("ACGTRYKMSWBDHVN") == "NBDHVWSKMRYACGT"


def test_orient_sequence_returns_forward_sequence_unchanged() -> None:
    assert orient_sequence("ACGTN", "forward") == "ACGTN"


def test_orient_sequence_reverse_complements_reverse_complement_records() -> None:
    assert orient_sequence("ACGTN", "reverse_complement") == "NACGT"


def test_orient_qualities_reverses_reverse_complement_records() -> None:
    assert orient_qualities([10, 20, 30, 40], "reverse_complement") == [40, 30, 20, 10]


def test_materialize_oriented_record_applies_orientation_to_sequence_and_qualities() -> (
    None
):
    record = make_record(
        sequence="AAGTC",
        qualities=[30, 31, 32, 33, 34],
        orientation="reverse_complement",
    )

    oriented_record = materialize_oriented_record(record)

    assert oriented_record.sequence == "GACTT"
    assert oriented_record.qualities == [34, 33, 32, 31, 30]
    assert oriented_record.orientation == "reverse_complement"


def test_orient_left_right_values_swaps_reverse_complement_display_values() -> None:
    assert orient_left_right_values(2, 5, "forward") == (2, 5)
    assert orient_left_right_values(2, 5, "reverse_complement") == (5, 2)


def test_orient_trim_config_for_display_swaps_left_and_right_for_reverse_complement() -> (
    None
):
    display_trim_config = orient_trim_config_for_display(
        TrimConfig(
            left_trim=2,
            right_trim=5,
            min_length=25,
            quality_trim_enabled=True,
            error_probability_cutoff=0.001,
        ),
        "reverse_complement",
    )

    assert display_trim_config == TrimConfig(
        left_trim=5,
        right_trim=2,
        min_length=25,
        quality_trim_enabled=True,
        error_probability_cutoff=0.001,
    )


def test_raw_trim_config_from_display_restores_raw_storage_orientation() -> None:
    raw_trim_config = raw_trim_config_from_display(
        TrimConfig(
            left_trim=5,
            right_trim=2,
            min_length=25,
            quality_trim_enabled=True,
            error_probability_cutoff=0.001,
        ),
        "reverse_complement",
    )

    assert raw_trim_config == TrimConfig(
        left_trim=2,
        right_trim=5,
        min_length=25,
        quality_trim_enabled=True,
        error_probability_cutoff=0.001,
    )
