from __future__ import annotations

from typing import Any

import pytest

from abi_sauce.models import SequenceUpload
from abi_sauce.parsers import abi


class FakeBioRecord:
    def __init__(self) -> None:
        self.id = "trace_001"
        self.name = "trace_001"
        self.description = "synthetic ABI record"
        self.seq = "ACGTN"
        self.annotations = {
            "machine_model": "SeqStudio",
            "sample_well": "A01",
            "run_start": "2026-03-11T10:00:00",
            "run_finish": "2026-03-11T10:30:00",
            "abif_raw": {
                "DATA9": [10, 20, 30],
                "DATA10": [5, 15, 25],
                "DATA11": [1, 2, 3],
                "DATA12": [7, 8, 9],
                "PLOC2": [100, 200, 300, 400, 500],
                "FWO_1": b"GATC\x00",
            },
        }
        self.letter_annotations = {
            "phred_quality": [40, 39, 38, 37, 10],
        }


def test_parse_ab1_upload_rejects_wrong_suffix() -> None:
    upload = SequenceUpload(filename="not_abi.fasta", content=b">x\nACGT\n")

    with pytest.raises(ValueError, match=r"Expected an \.ab1 or \.abi file"):
        abi.parse_ab1_upload(upload)


def test_parse_ab1_upload_returns_normalized_sequence_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_record = FakeBioRecord()

    def fake_read(*args: Any, **kwargs: Any) -> FakeBioRecord:
        return fake_record

    monkeypatch.setattr(abi.SeqIO, "read", fake_read)

    upload = SequenceUpload(filename="trace.ab1", content=b"fake-ab1-bytes")
    record = abi.parse_ab1_upload(upload)

    assert record.record_id == "trace_001"
    assert record.name == "trace_001"
    assert record.sequence == "ACGTN"
    assert record.source_format == "abi"

    assert record.qualities == [40, 39, 38, 37, 10]
    assert record.trace_data is not None
    assert sorted(record.trace_data.channels) == ["DATA10", "DATA11", "DATA12", "DATA9"]
    assert record.trace_data.base_positions == [100, 200, 300, 400, 500]
    assert record.trace_data.channel_order == "GATC"

    assert record.annotations["source_filename"] == "trace.ab1"
    assert record.annotations["machine_model"] == "SeqStudio"


def test_extract_trace_data_returns_none_when_no_channels_present() -> None:
    trace_data = abi._extract_trace_data({})
    assert trace_data is None
