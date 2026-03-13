from __future__ import annotations

from abi_sauce.chromatogram import build_chromatogram_view
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.trimming import TrimConfig, TrimResult, trim_sequence_record


def make_record(
    *,
    sequence: str = "ACGT",
    qualities: list[int] | None = None,
    trace_data: TraceData | None = None,
    name: str = "trace_001",
) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic trace",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=trace_data,
    )


def make_trace_data(
    *,
    channel_order: str | None = "GATC",
    channel_lengths: dict[str, int] | None = None,
    base_positions: list[int] | None = None,
) -> TraceData:
    channel_lengths = channel_lengths or {
        "DATA9": 50,
        "DATA10": 50,
        "DATA11": 50,
        "DATA12": 50,
    }
    channels = {
        data_key: list(range(length)) for data_key, length in channel_lengths.items()
    }
    return TraceData(
        channels=channels,
        base_positions=base_positions or [1, 3, 5, 7],
        channel_order=channel_order,
    )


def test_build_chromatogram_view_returns_not_renderable_without_trace_data() -> None:
    view = build_chromatogram_view(make_record(trace_data=None))

    assert view.is_renderable is False
    assert view.render_failure_reason == "missing_trace_data"
    assert view.channels == ()
    assert view.base_calls == ()
    assert view.quality_segments == ()


def test_build_chromatogram_view_renders_traces_and_base_calls_without_qualities() -> (
    None
):
    view = build_chromatogram_view(
        make_record(
            sequence="ACGT",
            qualities=None,
            trace_data=make_trace_data(),
        )
    )

    assert view.is_renderable is True
    assert view.trace_length == 50
    assert view.x_values == tuple(range(50))
    assert tuple(channel.base for channel in view.channels) == ("G", "A", "T", "C")
    assert tuple(base_call.base for base_call in view.base_calls) == tuple("ACGT")
    assert tuple(base_call.position for base_call in view.base_calls) == (1, 3, 5, 7)
    assert view.quality_segments == ()
    assert view.has_quality_overlay is False


def test_build_chromatogram_view_builds_quality_segments_from_peak_midpoints() -> None:
    view = build_chromatogram_view(
        make_record(
            sequence="ACGT",
            qualities=[40, 32, 18, 9],
            trace_data=make_trace_data(base_positions=[2, 4, 6, 8]),
        )
    )

    assert tuple(segment.left for segment in view.quality_segments) == (
        1.0,
        3.0,
        5.0,
        7.0,
    )
    assert tuple(segment.right for segment in view.quality_segments) == (
        3.0,
        5.0,
        7.0,
        9.0,
    )
    assert tuple(segment.center for segment in view.quality_segments) == (
        2.0,
        4.0,
        6.0,
        8.0,
    )
    assert tuple(segment.width for segment in view.quality_segments) == (
        2.0,
        2.0,
        2.0,
        2.0,
    )
    assert tuple(segment.quality for segment in view.quality_segments) == (
        40,
        32,
        18,
        9,
    )
    assert view.has_quality_overlay is True


def test_build_chromatogram_view_computes_trim_boundaries_from_base_midpoints() -> None:
    record = make_record(
        sequence="ACGT",
        qualities=[40, 41, 42, 43],
        trace_data=make_trace_data(base_positions=[10, 20, 30, 40]),
    )
    trim_result = TrimResult(
        record=record,
        original_length=4,
        trimmed_length=2,
        bases_removed_left=1,
        bases_removed_right=1,
        passed_min_length=True,
    )

    view = build_chromatogram_view(record, trim_result)

    assert view.trim_boundaries.left == 15.0
    assert view.trim_boundaries.right == 35.0


def test_build_chromatogram_view_falls_back_when_channel_order_is_missing_or_invalid() -> (
    None
):
    missing_order_view = build_chromatogram_view(
        make_record(trace_data=make_trace_data(channel_order=None))
    )
    invalid_order_view = build_chromatogram_view(
        make_record(trace_data=make_trace_data(channel_order="XYZQ"))
    )

    assert tuple(channel.base for channel in missing_order_view.channels) == (
        "G",
        "A",
        "T",
        "C",
    )
    assert tuple(channel.base for channel in invalid_order_view.channels) == (
        "G",
        "A",
        "T",
        "C",
    )


def test_build_chromatogram_view_sanitizes_channel_lengths_and_out_of_range_positions() -> (
    None
):
    view = build_chromatogram_view(
        make_record(
            sequence="ACGT",
            qualities=[40, 41, 42, 43],
            trace_data=make_trace_data(
                channel_lengths={
                    "DATA9": 9,
                    "DATA10": 7,
                    "DATA11": 8,
                    "DATA12": 6,
                },
                base_positions=[1, 3, 5, 8],
            ),
        )
    )

    assert view.trace_length == 6
    assert tuple(len(channel.signal) for channel in view.channels) == (6, 6, 6, 6)
    assert tuple(base_call.position for base_call in view.base_calls) == (1, 3, 5)
    assert tuple(segment.left for segment in view.quality_segments) == (0.0, 2.0, 4.0)
    assert tuple(segment.right for segment in view.quality_segments) == (2.0, 4.0, 5.0)


def test_build_chromatogram_view_keeps_renderable_trace_when_trim_fails_min_length() -> (
    None
):
    record = make_record(
        sequence="ACGT",
        qualities=[40, 39, 38, 37],
        trace_data=make_trace_data(base_positions=[2, 4, 6, 8]),
    )
    trim_result = trim_sequence_record(
        record,
        TrimConfig(left_trim=1, right_trim=1, min_length=5),
    )

    view = build_chromatogram_view(record, trim_result)

    assert trim_result.passed_min_length is False
    assert view.is_renderable is True
    assert view.trim_boundaries.left == 3.0
    assert view.trim_boundaries.right == 7.0
