from __future__ import annotations

from dataclasses import replace

from abi_sauce.models import SequenceOrientation, SequenceRecord
from abi_sauce.trimming import TrimConfig

_IUPAC_COMPLEMENT_TRANSLATION = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def complement_sequence(sequence: str) -> str:
    """Return the nucleotide complement of a sequence without reversing it."""
    return sequence.translate(_IUPAC_COMPLEMENT_TRANSLATION)


def complement_base(base: str) -> str:
    """Return the complemented base symbol for one base-call label."""
    return complement_sequence(base)


def reverse_complement_sequence(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return complement_sequence(sequence)[::-1]


def orient_sequence(sequence: str, orientation: SequenceOrientation) -> str:
    """Return a sequence in the requested display/export orientation."""
    if orientation == "forward":
        return sequence
    return reverse_complement_sequence(sequence)


def orient_qualities(
    qualities: list[int] | None,
    orientation: SequenceOrientation,
) -> list[int] | None:
    """Return per-base qualities in the requested display/export orientation."""
    if qualities is None:
        return None
    if orientation == "forward":
        return qualities
    return list(reversed(qualities))


def orient_left_right_values(
    left_value: int,
    right_value: int,
    orientation: SequenceOrientation,
) -> tuple[int, int]:
    """Return one left/right value pair in display-relative orientation."""
    if orientation == "forward":
        return (left_value, right_value)
    return (right_value, left_value)


def orient_trim_config_for_display(
    config: TrimConfig,
    orientation: SequenceOrientation,
) -> TrimConfig:
    """Return one trim config with left/right trims expressed in display space."""
    display_left_trim, display_right_trim = orient_left_right_values(
        config.left_trim,
        config.right_trim,
        orientation,
    )
    return replace(
        config,
        left_trim=display_left_trim,
        right_trim=display_right_trim,
    )


def raw_trim_config_from_display(
    config: TrimConfig,
    orientation: SequenceOrientation,
) -> TrimConfig:
    """Convert one display-relative trim config back into raw storage space."""
    raw_left_trim, raw_right_trim = orient_left_right_values(
        config.left_trim,
        config.right_trim,
        orientation,
    )
    return replace(
        config,
        left_trim=raw_left_trim,
        right_trim=raw_right_trim,
    )


def materialize_oriented_record(record: SequenceRecord) -> SequenceRecord:
    """Return a sequence record with orientation applied to sequence/qualities."""
    if record.orientation == "forward":
        return record

    return replace(
        record,
        sequence=orient_sequence(record.sequence, record.orientation),
        qualities=orient_qualities(record.qualities, record.orientation),
    )
