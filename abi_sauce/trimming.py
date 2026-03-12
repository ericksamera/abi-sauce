from __future__ import annotations

from dataclasses import dataclass, replace

from abi_sauce.models import SequenceRecord


@dataclass(frozen=True, slots=True)
class TrimConfig:
    """Configuration for simple sequence trimming."""

    left_trim: int = 0
    right_trim: int = 0
    min_length: int = 1


@dataclass(frozen=True, slots=True)
class TrimResult:
    """Result of trimming a sequence record."""

    record: SequenceRecord
    original_length: int
    trimmed_length: int
    bases_removed_left: int
    bases_removed_right: int
    passed_min_length: bool


def trim_sequence_record(
    record: SequenceRecord,
    config: TrimConfig,
) -> TrimResult:
    """Trim a sequence record using fixed left/right clipping."""
    if config.left_trim < 0:
        raise ValueError("left_trim must be >= 0")
    if config.right_trim < 0:
        raise ValueError("right_trim must be >= 0")
    if config.min_length < 0:
        raise ValueError("min_length must be >= 0")

    original_sequence = record.sequence
    original_length = len(original_sequence)

    start = min(config.left_trim, original_length)
    end = original_length - config.right_trim

    if end < start:
        end = start

    trimmed_sequence = original_sequence[start:end]

    trimmed_qualities: list[int] | None = None
    if record.qualities is not None:
        trimmed_qualities = record.qualities[start:end]

    trimmed_record = replace(
        record,
        sequence=trimmed_sequence,
        qualities=trimmed_qualities,
    )

    trimmed_length = len(trimmed_sequence)

    return TrimResult(
        record=trimmed_record,
        original_length=original_length,
        trimmed_length=trimmed_length,
        bases_removed_left=start,
        bases_removed_right=original_length - end,
        passed_min_length=trimmed_length >= config.min_length,
    )
