from __future__ import annotations

from typing import Any

import pytest

from abi_sauce.models import SequenceUpload
from abi_sauce.parsers import abi


class FakeBioRecord:
    def __init__(
        self,
        *,
        annotations: dict[str, Any] | None = None,
        letter_annotations: dict[str, Any] | None = None,
    ) -> None:
        self.id = "trace_001"
        self.name = "trace_001"
        self.description = "synthetic ABI record"
        self.seq = "ACGTN"
        self.annotations = annotations or {}
        self.letter_annotations = letter_annotations or {}


def make_fake_bio_record(
    *,
    include_abif_raw: bool = True,
    include_phred_quality: bool = True,
    use_ploc2: bool = True,
    use_ploc1: bool = False,
    include_fwo_1: bool = True,
    include_qc_scores: bool = True,
) -> FakeBioRecord:
    annotations: dict[str, Any] = {
        "machine_model": "SeqStudio",
        "sample_well": "A01",
        "run_start": "2026-03-11T10:00:00",
        "run_finish": "2026-03-11T10:30:00",
    }

    if include_abif_raw:
        abif_raw: dict[str, Any] = {
            "DATA9": [10, 20, 30],
            "DATA10": [5, 15, 25],
            "DATA11": [1, 2, 3],
            "DATA12": [7, 8, 9],
        }
        if use_ploc2:
            abif_raw["PLOC2"] = [100, 200, 300, 400, 500]
        if use_ploc1:
            abif_raw["PLOC1"] = [101, 201, 301, 401, 501]
        if include_fwo_1:
            abif_raw["FWO_1"] = b"GATC\x00"
        if include_qc_scores:
            abif_raw["TrSc1"] = 31
            abif_raw["PuSc1"] = 23
            abif_raw["CRLn1"] = 456

        annotations["abif_raw"] = abif_raw

    letter_annotations: dict[str, Any] = {}
    if include_phred_quality:
        letter_annotations["phred_quality"] = [40, 39, 38, 37, 10]

    return FakeBioRecord(
        annotations=annotations,
        letter_annotations=letter_annotations,
    )


def patch_seqio_read(
    monkeypatch: pytest.MonkeyPatch,
    *,
    record: FakeBioRecord | None = None,
    side_effect: Exception | None = None,
) -> None:
    def fake_read(*args: Any, **kwargs: Any) -> FakeBioRecord:
        if side_effect is not None:
            raise side_effect
        assert record is not None
        return record

    monkeypatch.setattr(abi.SeqIO, "read", fake_read)


def test_parse_ab1_upload_rejects_wrong_suffix() -> None:
    upload = SequenceUpload(filename="not_abi.fasta", content=b">x\nACGT\n")

    with pytest.raises(ValueError, match=r"Expected an \.ab1 or \.abi file"):
        abi.parse_ab1_upload(upload)


def test_parse_ab1_upload_returns_normalized_sequence_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(monkeypatch, record=make_fake_bio_record())

    upload = SequenceUpload(filename="trace.ab1", content=b"fake-ab1-bytes")
    record = abi.parse_ab1_upload(upload)

    assert record.record_id == "trace_001"
    assert record.name == "trace_001"
    assert record.sequence == "ACGTN"
    assert record.source_format == "abi"
    assert record.orientation == "forward"

    assert record.qualities == [40, 39, 38, 37, 10]
    assert record.trace_data is not None
    assert sorted(record.trace_data.channels) == ["DATA10", "DATA11", "DATA12", "DATA9"]
    assert record.trace_data.base_positions == [100, 200, 300, 400, 500]
    assert record.trace_data.channel_order == "GATC"

    assert record.annotations["source_filename"] == "trace.ab1"
    assert record.annotations["machine_model"] == "SeqStudio"
    assert record.annotations["trace_score"] == 31
    assert record.annotations["pup_score"] == 23
    assert record.annotations["crl_score"] == 456


def test_extract_trace_data_returns_none_when_no_channels_present() -> None:
    trace_data = abi._extract_trace_data({})
    assert trace_data is None


def test_parse_ab1_upload_accepts_uppercase_and_mixed_case_suffixes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(monkeypatch, record=make_fake_bio_record())

    upper = SequenceUpload(filename="trace.AB1", content=b"fake")
    mixed = SequenceUpload(filename="trace.Abi", content=b"fake")

    upper_record = abi.parse_ab1_upload(upper)
    mixed_record = abi.parse_ab1_upload(mixed)

    assert upper_record.source_format == "abi"
    assert mixed_record.source_format == "abi"


def test_parse_ab1_upload_succeeds_when_phred_quality_is_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(
        monkeypatch,
        record=make_fake_bio_record(include_phred_quality=False),
    )

    upload = SequenceUpload(filename="trace.ab1", content=b"fake")
    record = abi.parse_ab1_upload(upload)

    assert record.sequence == "ACGTN"
    assert record.qualities is None


def test_parse_ab1_upload_succeeds_when_abif_raw_is_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(
        monkeypatch,
        record=make_fake_bio_record(include_abif_raw=False),
    )

    upload = SequenceUpload(filename="trace.ab1", content=b"fake")
    record = abi.parse_ab1_upload(upload)

    assert record.sequence == "ACGTN"
    assert record.trace_data is None


def test_extract_trace_data_falls_back_to_ploc1_when_ploc2_is_missing() -> None:
    trace_data = abi._extract_trace_data(
        {
            "DATA9": [10, 20, 30],
            "DATA10": [5, 15, 25],
            "DATA11": [1, 2, 3],
            "DATA12": [7, 8, 9],
            "PLOC1": [101, 201, 301, 401, 501],
            "FWO_1": b"GATC\x00",
        }
    )

    assert trace_data is not None
    assert trace_data.base_positions == [101, 201, 301, 401, 501]


def test_extract_trace_data_succeeds_when_fwo_1_is_missing() -> None:
    trace_data = abi._extract_trace_data(
        {
            "DATA9": [10, 20, 30],
            "DATA10": [5, 15, 25],
            "DATA11": [1, 2, 3],
            "DATA12": [7, 8, 9],
            "PLOC2": [100, 200, 300, 400, 500],
        }
    )

    assert trace_data is not None
    assert trace_data.channel_order is None


def test_parse_ab1_upload_sets_missing_qc_metrics_to_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(
        monkeypatch,
        record=make_fake_bio_record(include_qc_scores=False),
    )

    upload = SequenceUpload(filename="trace.ab1", content=b"fake")
    record = abi.parse_ab1_upload(upload)

    assert record.annotations["trace_score"] is None
    assert record.annotations["pup_score"] is None
    assert record.annotations["crl_score"] is None


def test_parse_ab1_upload_wraps_low_level_parse_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    patch_seqio_read(monkeypatch, side_effect=ValueError("corrupt ABI bytes"))

    upload = SequenceUpload(filename="broken.ab1", content=b"not-a-real-ab1")

    with pytest.raises(
        abi.AbiParseError,
        match=r"Failed to parse ABI file: broken\.ab1",
    ):
        abi.parse_ab1_upload(upload)
