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
    error_probability_cutoff: float = 0.01


@dataclass(frozen=True, slots=True)
class TrimResult:
    """Result of trimming a sequence record."""

    record: SequenceRecord
    original_length: int
    trimmed_length: int
    bases_removed_left: int
    bases_removed_right: int
    passed_min_length: bool
    fixed_bases_removed_left: int = 0
    fixed_bases_removed_right: int = 0
    quality_bases_removed_left: int = 0
    quality_bases_removed_right: int = 0


def trim_sequence_record(
    record: SequenceRecord,
    config: TrimConfig,
) -> TrimResult:
    """Trim a sequence record using fixed and optional Mott quality trimming."""
    _validate_trim_config(config)
    _validate_record_for_trimming(record)

    original_sequence = record.sequence
    original_length = len(original_sequence)

    quality_start = 0
    quality_end = original_length

    if config.quality_trim_enabled and record.qualities is not None:
        quality_start, quality_end = _mott_trim_bounds(
            record.qualities,
            config.error_probability_cutoff,
        )

    available_after_quality = max(quality_end - quality_start, 0)
    fixed_bases_removed_left = min(config.left_trim, available_after_quality)
    fixed_bases_removed_right = min(
        config.right_trim,
        max(available_after_quality - fixed_bases_removed_left, 0),
    )

    start = quality_start + fixed_bases_removed_left
    end = quality_end - fixed_bases_removed_right

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
        fixed_bases_removed_left=fixed_bases_removed_left,
        fixed_bases_removed_right=fixed_bases_removed_right,
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
    if not 0 <= config.error_probability_cutoff <= 1:
        raise ValueError("error_probability_cutoff must be between 0 and 1")


def _validate_record_for_trimming(record: SequenceRecord) -> None:
    """Validate record fields required for trimming."""
    if record.qualities is not None and len(record.qualities) != len(record.sequence):
        raise ValueError(
            "record.qualities must have the same length as record.sequence"
        )


def _mott_trim_bounds(
    qualities: list[int],
    cutoff: float,
) -> tuple[int, int]:
    """Return trimmed ``[start, end)`` bounds using modified Mott trimming.

    Each base contributes ``cutoff - P(error)`` where
    ``P(error) = 10 ** (-Q / 10)`` for PHRED score ``Q``. The retained region is
    the maximum-scoring contiguous segment, with negative cumulative scores reset
    to zero.

    This matches the R variant that does not forcibly discard the first base.
    """
    if not qualities:
        return 0, 0

    best_start = 0
    best_end = 0
    best_score = 0.0

    current_start = 0
    current_score = 0.0

    for index, quality in enumerate(qualities):
        base_score = cutoff - _phred_error_probability(quality)
        current_score += base_score

        if current_score <= 0:
            current_score = 0.0
            current_start = index + 1
            continue

        if current_score > best_score:
            best_score = current_score
            best_start = current_start
            best_end = index + 1

    return best_start, best_end


def _phred_error_probability(quality: int) -> float:
    """Convert a PHRED score into an estimated base-call error probability."""
    return 10 ** (-quality / 10.0)
