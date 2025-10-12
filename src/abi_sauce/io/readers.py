from __future__ import annotations
from io import StringIO, BytesIO
from typing import List, Dict, Any
from uuid import uuid4

from Bio import SeqIO

from abi_sauce.models import SequenceAsset, TraceAsset


def _mk_id() -> str:
    return uuid4().hex


def _decode_text(raw: bytes) -> str:
    """Decode bytes to text for FASTA/GenBank/ApE. UTF-8 first, fallback to Latin-1."""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


# ---------- FASTA ----------
def read_fasta(name: str, raw: bytes, size: int, checksum: str) -> List[SequenceAsset]:
    # Biopython expects text-mode handle for text formats like FASTA
    handle = StringIO(_decode_text(raw))
    assets: List[SequenceAsset] = []
    for rec in SeqIO.parse(handle, "fasta"):
        assets.append(
            SequenceAsset(
                id=_mk_id(),
                name=name,
                ext="fasta",
                size=size,
                checksum_md5=checksum,
                sequence=str(rec.seq),
                description=rec.description or rec.id,
                meta={"id": rec.id},
            )
        )
    return assets


# ---------- GenBank / ApE ----------
def read_genbank_like(
    name: str, raw: bytes, size: int, checksum: str, ext_hint: str
) -> List[SequenceAsset]:
    """Read GenBank (and ApE, which is GenBank-like) as sequences with features."""
    # Text-mode handle required for 'genbank'
    handle = StringIO(_decode_text(raw))
    assets: List[SequenceAsset] = []
    for rec in SeqIO.parse(handle, "genbank"):
        feats = []
        for f in rec.features:
            feats.append(
                {
                    "type": f.type,
                    "location": str(f.location),
                    "qualifiers": {k: v for k, v in (f.qualifiers or {}).items()},
                }
            )
        assets.append(
            SequenceAsset(
                id=_mk_id(),
                name=name,
                ext=ext_hint,
                size=size,
                checksum_md5=checksum,
                sequence=str(rec.seq),
                description=rec.description or rec.id,
                meta={"id": rec.id, "annotations": dict(rec.annotations or {})},
                features=feats,
            )
        )
    return assets


# ---------- AB1 ----------
def read_ab1(name: str, raw: bytes, size: int, checksum: str) -> List[TraceAsset]:
    # Binary-mode handle for AB1 (ABIF)
    handle = BytesIO(raw)
    rec = SeqIO.read(handle, "abi")  # raises on failure

    abif_raw: Dict[str, Any] = dict(rec.annotations.get("abif_raw", {}))

    # These keys are common across many instruments.
    # Raw channel arrays are often DATA9..DATA12, base order given by FWO_ (e.g., b"GATC").
    # PLOC2 gives base positions (peak indices), PBAS2 gives called sequence, PCON2 gives quality.
    fwo = abif_raw.get("FWO_", b"GATC")
    order = (
        fwo.decode(errors="ignore") if isinstance(fwo, (bytes, bytearray)) else str(fwo)
    )
    data_keys = ["DATA9", "DATA10", "DATA11", "DATA12"]

    channels = {b: None for b in "ACGT"}
    for base, key in zip(order, data_keys):
        arr = abif_raw.get(key)
        if arr is not None:
            try:
                channels[base] = list(map(int, arr))
            except Exception:
                channels[base] = None

    seq = str(rec.seq) if getattr(rec, "seq", None) else None

    # Qualities (PHRED)
    quals = None
    if hasattr(rec, "letter_annotations") and "phred_quality" in rec.letter_annotations:
        quals = list(map(int, rec.letter_annotations["phred_quality"]))
    elif "PCON2" in abif_raw:
        try:
            quals = list(map(int, abif_raw["PCON2"]))
        except Exception:
            pass

    # Base positions (peaks)
    ploc = None
    if "PLOC2" in abif_raw:
        try:
            ploc = list(map(int, abif_raw["PLOC2"]))
        except Exception:
            pass

    asset = TraceAsset(
        id=_mk_id(),
        name=name,
        ext="ab1",
        size=size,
        checksum_md5=checksum,
        sequence=seq or None,
        qualities=quals,
        base_positions=ploc,
        channels=channels,
        meta={
            "abif_keys": list(abif_raw.keys())[:40]
        },  # keep lightweight; avoid dumping huge dict
    )
    return [asset]
