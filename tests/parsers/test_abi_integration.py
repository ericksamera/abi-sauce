from __future__ import annotations

from abi_sauce.parsers.abi import parse_ab1_upload


def test_parse_real_ab1_file(real_ab1_upload) -> None:
    record = parse_ab1_upload(real_ab1_upload)

    assert record.source_format == "abi"
    assert record.sequence
    assert len(record.sequence) > 0

    if record.qualities is not None:
        assert len(record.qualities) == len(record.sequence)

    assert record.trace_data is not None
    assert len(record.trace_data.channels) == 4

    for signal in record.trace_data.channels.values():
        assert len(signal) > 0

    assert len(record.trace_data.base_positions) > 0
    assert record.trace_data.base_positions == sorted(record.trace_data.base_positions)
    assert record.trace_data.base_positions[0] >= 0
    assert record.trace_data.base_positions[-1] > record.trace_data.base_positions[0]
