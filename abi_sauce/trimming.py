from __future__ import annotations

from dataclasses import dataclass, replace

from abi_sauce.models import SequenceRecord


@dataclass(frozen=True, slots=True)
class TrimConfig:
    """Configuration for fixed and quality-aware sequence trimming."""

    left_trim: int = 0
    right_trim: int = 0
    min_length: int = 1
    quality_trim_enabled: bool = False
    quality_threshold: int = 20


@dataclass(frozen=True, slots=True)
class TrimResult:
    """Result of trimming a sequence record."""

    record: SequenceRecord
    original_length: int
    trimmed_length: int
    bases_removed_left: int
    bases_removed_right: int
    passed_min_length: bool
    quality_bases_removed_left: int = 0
    quality_bases_removed_right: int = 0


def trim_sequence_record(
    record: SequenceRecord,
    config: TrimConfig,
) -> TrimResult:
    """Trim a sequence record using fixed and optional quality-aware clipping."""
    _validate_trim_config(config)

    original_sequence = record.sequence
    original_length = len(original_sequence)

    quality_start = 0
    quality_end = original_length

    if config.quality_trim_enabled and record.qualities is not None:
        quality_start, quality_end = _quality_trim_bounds(
            record.qualities,
            config.quality_threshold,
        )

    start = min(quality_start + config.left_trim, original_length)
    end = quality_end - config.right_trim

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
        quality_bases_removed_left=quality_start,
        quality_bases_removed_right=original_length - quality_end,
    )


def _validate_trim_config(config: TrimConfig) -> None:
    """Validate trimming configuration values."""
    if config.left_trim < 0:
        raise ValueError("left_trim must be >= 0")
    if config.right_trim < 0:
        raise ValueError("right_trim must be >= 0")
    if config.min_length < 0:
        raise ValueError("min_length must be >= 0")
    if config.quality_threshold < 0:
        raise ValueError("quality_threshold must be >= 0")


def _quality_trim_bounds(
    qualities: list[int],
    threshold: int,
) -> tuple[int, int]:
    """Return trimmed [start, end) bounds after end-only quality clipping."""
    start = 0
    end = len(qualities)

    while start < end and qualities[start] < threshold:
        start += 1

    while end > start and qualities[end - 1] < threshold:
        end -= 1

    return start, end
