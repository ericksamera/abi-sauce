from __future__ import annotations

from io import BytesIO
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord as BioSeqRecord

from abi_sauce.exceptions import AbiParseError
from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData

_TRACE_KEYS = ("DATA9", "DATA10", "DATA11", "DATA12")


def parse_ab1_upload(upload: SequenceUpload) -> SequenceRecord:
    """Parse an ABI/AB1 upload into a normalized SequenceRecord."""
    if upload.suffix not in {"ab1", "abi"}:
        raise ValueError(f"Expected an .ab1 or .abi file, got: {upload.filename}")

    try:
        bio_record = SeqIO.read(BytesIO(upload.content), "abi")
    except Exception as exc:
        raise AbiParseError(f"Failed to parse ABI file: {upload.filename}") from exc

    return _to_sequence_record(bio_record=bio_record, source_filename=upload.filename)


def _to_sequence_record(
    *,
    bio_record: BioSeqRecord,
    source_filename: str,
) -> SequenceRecord:
    raw_abif = bio_record.annotations.get("abif_raw")
    raw_annotations: dict[str, Any]
    if isinstance(raw_abif, dict):
        raw_annotations = raw_abif
    else:
        raw_annotations = {}
    qualities = _extract_qualities(bio_record)
    trace_data = _extract_trace_data(raw_annotations)

    return SequenceRecord(
        record_id=bio_record.id or source_filename,
        name=bio_record.name or bio_record.id or source_filename,
        description=bio_record.description,
        sequence=str(bio_record.seq),
        source_format="abi",
        qualities=qualities,
        trace_data=trace_data,
        annotations={
            "source_filename": source_filename,
            "machine_model": bio_record.annotations.get("machine_model"),
            "sample_well": bio_record.annotations.get("sample_well"),
            "run_start": bio_record.annotations.get("run_start"),
            "run_finish": bio_record.annotations.get("run_finish"),
        },
    )


def _extract_qualities(bio_record: BioSeqRecord) -> list[int] | None:
    """Extract PHRED qualities from a Biopython SeqRecord."""
    raw_qualities = bio_record.letter_annotations.get("phred_quality")
    if raw_qualities is None:
        return None
    return [int(value) for value in raw_qualities]


def _extract_trace_data(raw_annotations: dict[str, Any]) -> TraceData | None:
    """Extract ABI trace channels and peak locations from ABIF raw annotations."""
    channels = {
        key: [int(value) for value in raw_annotations.get(key, [])]
        for key in _TRACE_KEYS
        if key in raw_annotations
    }
    if not channels:
        return None

    base_positions_raw = (
        raw_annotations.get("PLOC2") or raw_annotations.get("PLOC1") or []
    )
    channel_order = _decode_text(raw_annotations.get("FWO_1"))

    return TraceData(
        channels=channels,
        base_positions=[int(value) for value in base_positions_raw],
        channel_order=channel_order,
    )


def _decode_text(value: Any) -> str | None:
    """Convert an ABIF text-like field into a Python string."""
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").strip("\x00")
    return str(value)
