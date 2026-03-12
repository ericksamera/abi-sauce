from __future__ import annotations

import pytest

from abi_sauce.models import SequenceRecord
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    sequence: str = "ACGTACGT",
    qualities: list[int] | None = None,
) -> SequenceRecord:
    return SequenceRecord(
        record_id="r1",
        name="r1",
        description="test",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
    )


def test_trim_sequence_record_noop() -> None:
    record = make_record(sequence="ACGT")
    result = trim_sequence_record(record, TrimConfig())

    assert result.record.sequence == "ACGT"
    assert result.original_length == 4
    assert result.trimmed_length == 4
    assert result.bases_removed_left == 0
    assert result.bases_removed_right == 0
    assert result.passed_min_length is True


def test_trim_sequence_record_left_trim() -> None:
    record = make_record(sequence="ACGTACGT")
    result = trim_sequence_record(record, TrimConfig(left_trim=2))

    assert result.record.sequence == "GTACGT"
    assert result.bases_removed_left == 2
    assert result.bases_removed_right == 0


def test_trim_sequence_record_right_trim() -> None:
    record = make_record(sequence="ACGTACGT")
    result = trim_sequence_record(record, TrimConfig(right_trim=3))

    assert result.record.sequence == "ACGTA"
    assert result.bases_removed_left == 0
    assert result.bases_removed_right == 3


def test_trim_sequence_record_left_and_right_trim() -> None:
    record = make_record(sequence="ACGTACGT")
    result = trim_sequence_record(
        record,
        TrimConfig(left_trim=2, right_trim=2),
    )

    assert result.record.sequence == "GTAC"
    assert result.trimmed_length == 4


def test_trim_sequence_record_trims_qualities_in_sync() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[10, 11, 12, 13, 14, 15, 16, 17],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(left_trim=1, right_trim=2),
    )

    assert result.record.sequence == "CGTAC"
    assert result.record.qualities == [11, 12, 13, 14, 15]


def test_trim_sequence_record_flags_min_length_failure() -> None:
    record = make_record(sequence="ACGT")
    result = trim_sequence_record(
        record,
        TrimConfig(left_trim=1, right_trim=1, min_length=5),
    )

    assert result.record.sequence == "CG"
    assert result.passed_min_length is False


def test_trim_sequence_record_handles_overtrim() -> None:
    record = make_record(sequence="ACGT")
    result = trim_sequence_record(
        record,
        TrimConfig(left_trim=3, right_trim=3),
    )

    assert result.record.sequence == ""
    assert result.trimmed_length == 0


def test_trim_sequence_record_rejects_negative_values() -> None:
    record = make_record()

    with pytest.raises(ValueError, match="left_trim must be >= 0"):
        trim_sequence_record(record, TrimConfig(left_trim=-1))


def test_trim_sequence_record_quality_trim_left_end() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[5, 8, 30, 35, 40, 41, 42, 43],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "GTACGT"
    assert result.record.qualities == [30, 35, 40, 41, 42, 43]
    assert result.quality_bases_removed_left == 2
    assert result.quality_bases_removed_right == 0


def test_trim_sequence_record_quality_trim_right_end() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[30, 31, 32, 33, 34, 10, 5, 3],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "ACGTA"
    assert result.record.qualities == [30, 31, 32, 33, 34]
    assert result.quality_bases_removed_left == 0
    assert result.quality_bases_removed_right == 3


def test_trim_sequence_record_quality_trim_both_ends() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[5, 8, 30, 35, 40, 25, 7, 3],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "GTAC"
    assert result.record.qualities == [30, 35, 40, 25]
    assert result.quality_bases_removed_left == 2
    assert result.quality_bases_removed_right == 2


def test_trim_sequence_record_quality_trim_keeps_low_quality_internal_bases() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[30, 30, 5, 30, 5, 30, 30, 30],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "ACGTACGT"
    assert result.record.qualities == [30, 30, 5, 30, 5, 30, 30, 30]
    assert result.quality_bases_removed_left == 0
    assert result.quality_bases_removed_right == 0


def test_trim_sequence_record_quality_trim_is_skipped_when_qualities_missing() -> None:
    record = make_record(sequence="ACGTACGT", qualities=None)
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "ACGTACGT"
    assert result.record.qualities is None
    assert result.quality_bases_removed_left == 0
    assert result.quality_bases_removed_right == 0


def test_trim_sequence_record_quality_trim_and_fixed_trim_combine() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[5, 8, 30, 35, 40, 25, 7, 3],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            left_trim=1,
            right_trim=1,
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == "TA"
    assert result.record.qualities == [35, 40]
    assert result.quality_bases_removed_left == 2
    assert result.quality_bases_removed_right == 2
    assert result.bases_removed_left == 3
    assert result.bases_removed_right == 3


def test_trim_sequence_record_quality_trim_all_bases_below_threshold() -> None:
    record = make_record(
        sequence="ACGT",
        qualities=[5, 6, 7, 8],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
        ),
    )

    assert result.record.sequence == ""
    assert result.record.qualities == []
    assert result.trimmed_length == 0
    assert result.quality_bases_removed_left == 4
    assert result.quality_bases_removed_right == 0


def test_trim_sequence_record_quality_trim_can_fail_min_length() -> None:
    record = make_record(
        sequence="ACGTACGT",
        qualities=[5, 8, 30, 35, 40, 25, 7, 3],
    )
    result = trim_sequence_record(
        record,
        TrimConfig(
            quality_trim_enabled=True,
            quality_threshold=20,
            min_length=5,
        ),
    )

    assert result.record.sequence == "GTAC"
    assert result.trimmed_length == 4
    assert result.passed_min_length is False


def test_trim_sequence_record_rejects_negative_quality_threshold() -> None:
    record = make_record()

    with pytest.raises(ValueError, match="quality_threshold must be >= 0"):
        trim_sequence_record(
            record,
            TrimConfig(
                quality_trim_enabled=True,
                quality_threshold=-1,
            ),
        )
