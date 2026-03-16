from __future__ import annotations

import pytest

from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData
from abi_sauce.services.batch import (
    ParsedBatch,
    apply_trim_configs,
    build_batch_signature,
)
from abi_sauce.services.reference_alignment import compute_reference_alignment
from abi_sauce.trimming import TrimConfig


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    trace_length = 100
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic record",
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
    uploads = (SequenceUpload(filename="trace.ab1", content=b"trace"),)
    parsed_batch = ParsedBatch(
        uploads=uploads,
        parsed_records={
            "trace.ab1": make_record(
                name="trace",
                sequence="AACCGGTT",
                qualities=[10, 20, 30, 40, 50, 40, 30, 20],
            )
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )
    return apply_trim_configs(
        parsed_batch,
        trim_configs_by_name={"trace.ab1": TrimConfig(left_trim=2, right_trim=2)},
    )


def test_compute_reference_alignment_returns_page_facing_alignment_state() -> None:
    prepared_batch = make_prepared_batch()

    computed_alignment = compute_reference_alignment(
        prepared_batch,
        source_filename="trace.ab1",
        reference_text=">ref\nCCGA\n",
        strand_policy="forward",
    )

    assert computed_alignment.source_filename == "trace.ab1"
    assert computed_alignment.alignment_result.reference_name == "ref"
    assert computed_alignment.alignment_result.mismatch_count == 1
    assert len(computed_alignment.event_rows) == 1
    assert computed_alignment.event_rows[0]["type"] == "mismatch"
    assert computed_alignment.chromatogram_view.is_renderable is True


def test_compute_reference_alignment_can_include_match_rows() -> None:
    prepared_batch = make_prepared_batch()

    computed_alignment = compute_reference_alignment(
        prepared_batch,
        source_filename="trace.ab1",
        reference_text=">ref\nCCGG\n",
        strand_policy="forward",
        include_matches=True,
    )

    assert len(computed_alignment.event_rows) == 4
    assert all(row["type"] == "match" for row in computed_alignment.event_rows)


def test_compute_reference_alignment_rejects_unknown_source_filename() -> None:
    prepared_batch = make_prepared_batch()

    with pytest.raises(KeyError, match="missing.ab1"):
        compute_reference_alignment(
            prepared_batch,
            source_filename="missing.ab1",
            reference_text=">ref\nCCGG\n",
        )
