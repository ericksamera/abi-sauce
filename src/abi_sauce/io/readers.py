from __future__ import annotations

from io import StringIO
from typing import Any, Protocol
from uuid import uuid4

from Bio import SeqIO

from abi_sauce.io.detect import FileKind
from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset
from abi_sauce.services.feature_ops import feature_from_biopython, features_to_dicts


def _mk_id() -> str:
    return uuid4().hex


def _decode_text(raw: bytes) -> str:
    """Decode bytes to text for FASTA/GenBank/ApE. UTF-8 first, fallback to Latin-1."""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


# ---------- FASTA ----------
def read_fasta(name: str, raw: bytes, size: int, checksum: str) -> list[SequenceAsset]:
    handle = StringIO(_decode_text(raw))
    assets: list[SequenceAsset] = []
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
) -> list[SequenceAsset]:
    """Read GenBank (and ApE, which is GenBank-like) as sequences with features."""
    handle = StringIO(_decode_text(raw))
    assets: list[SequenceAsset] = []
    for rec in SeqIO.parse(handle, "genbank"):
        # Convert to typed Feature, then to dicts (UI-friendly) for now.
        typed = [feature_from_biopython(f, len(rec.seq)) for f in rec.features]
        feats = features_to_dicts(typed)

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
def read_ab1(name: str, raw: bytes, size: int, checksum: str) -> list[TraceAsset]:
    from abi_sauce.services.ab1 import parse_ab1

    channels, seq, quals, ploc, meta = parse_ab1(raw)

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
        meta=meta,
    )
    return [asset]


# ---------- Registry ----------
class ReaderFn(Protocol):
    def __call__(
        self, name: str, raw: bytes, size: int, checksum: str, **kwargs: Any
    ) -> list[AssetBase]: ...


READERS: dict[FileKind, ReaderFn] = {
    "fasta": read_fasta,
    "ab1": read_ab1,
    # genbank & ape share the same parser with an ext hint
    "genbank": lambda name, raw, size, checksum, **_: read_genbank_like(
        name, raw, size, checksum, ext_hint="gbk"
    ),
    "ape": lambda name, raw, size, checksum, **_: read_genbank_like(
        name, raw, size, checksum, ext_hint="ape"
    ),
    # "unknown" intentionally omitted
}
