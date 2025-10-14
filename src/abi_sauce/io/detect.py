from __future__ import annotations

from typing import Literal, cast

FileKind = Literal["fasta", "genbank", "ape", "ab1", "unknown"]

FASTA_FIRST_CHARS = (b">",)
AB1_MAGIC = b"ABIF"  # ABIF header for ABI chromatogram files

EXT_MAP: dict[str, FileKind] = {
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbff": "genbank",
    ".ape": "ape",
    ".ab1": "ab1",
    ".abi": "ab1",
}


def by_extension(name: str) -> FileKind:
    lname = name.lower()
    for ext, kind in EXT_MAP.items():
        if lname.endswith(ext):
            return cast(FileKind, kind)
    return "unknown"


def sniff_bytes(name: str, data: bytes) -> FileKind:
    kind = by_extension(name)
    if kind == "ab1":
        return "ab1"
    if data.startswith(AB1_MAGIC):
        return "ab1"
    if data[:1] in FASTA_FIRST_CHARS:
        return "fasta"
    # fallback: many GenBank/ApE start with LOCUS/FEATURES/ORIGIN; keep ext-based mapping
    return kind


def identify(name: str, data: bytes) -> tuple[FileKind, str]:
    """Return (kind, normalized_ext)."""
    kind = sniff_bytes(name, data)
    # normalize ext for downstream: e.g., "ape" treated as genbank-like on read
    norm_ext = {
        "fasta": "fasta",
        "genbank": "gbk",
        "ape": "ape",
        "ab1": "ab1",
        "unknown": "",
    }[kind]
    return kind, norm_ext
