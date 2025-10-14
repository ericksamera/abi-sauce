# src/abi_sauce/services/sequence_ops.py
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

# Lightweight helpers; Biopython is optional here (import only when needed)

_COMP = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def revcomp(seq: str) -> str:
    """Reverse-complement DNA (keeps case, supports IUPAC)."""
    return (seq or "").translate(_COMP)[::-1]


def gc_content(seq: str) -> float:
    """GC fraction in [0,1] (ignores non-ACGT)."""
    s = (seq or "").upper()
    if not s:
        return 0.0
    gc = sum(1 for c in s if c in "GC")
    atgc = sum(1 for c in s if c in "ACGT")
    return (gc / atgc) if atgc else 0.0


def translate(
    seq: str, *, frame: int = 0, table: int = 11, to_stop: bool = False
) -> str:
    """Translate with Biopython if available; else a tiny fallback that won't support alt tables."""
    try:
        from Bio.Seq import Seq

        return str(Seq(seq[frame:]).translate(table=table, to_stop=to_stop))
    except Exception:
        CODON = {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "TAT": "Y",
            "TAC": "Y",
            "TAA": "*",
            "TAG": "*",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "TGT": "C",
            "TGC": "C",
            "TGA": "*",
            "TGG": "W",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }
        s = (seq or "").upper()[frame:]
        aa = []
        for i in range(0, len(s) - 2, 3):
            aa1 = CODON.get(s[i : i + 3], "X")
            if to_stop and aa1 == "*":
                break
            aa.append(aa1)
        return "".join(aa)


@dataclass(frozen=True)
class ORF:
    start: int  # 1-based inclusive on the forward strand
    end: int  # 1-based inclusive
    frame: int  # 0,1,2 for +; -1,-2,-3 for reverse
    strand: int  # +1 or -1
    aa_len: int
    prod: str


def find_orfs(seq: str, *, min_aa: int = 30, table: int = 11) -> list[ORF]:
    """Simple ORF finder on + and - strands; returns AA-products."""

    def _scan(s: str, base_frame: int, strand: int) -> list[ORF]:
        aa = translate(s, frame=0, table=table, to_stop=False)
        out: list[ORF] = []
        start = 0
        for i, aa_c in enumerate(aa + "*"):
            if aa_c == "M" and start == 0:
                start = i + 1  # 1-based aa
            if aa_c == "*" and start > 0:
                aa_len = i - (start - 1)
                if aa_len >= min_aa:
                    # Map aa coords back to nucleotide coords (1-based)
                    nt_start = (start - 1) * 3 + 1
                    nt_end = i * 3
                    if strand == +1:
                        out.append(
                            ORF(
                                nt_start + base_frame,
                                nt_end + base_frame,
                                base_frame,
                                +1,
                                aa_len,
                                aa[start - 1 : i],
                            )
                        )
                    else:
                        # on reverse we map relative to original sense
                        # we use forward coords by reflecting later
                        out.append(
                            ORF(
                                nt_start + base_frame,
                                nt_end + base_frame,
                                -(base_frame + 1),
                                -1,
                                aa_len,
                                aa[start - 1 : i],
                            )
                        )
                start = 0
        return out

    s = seq or ""
    outs: list[ORF] = []
    # forward
    for f in (0, 1, 2):
        outs.extend(_scan(s[f:], f, +1))
    # reverse
    rc = revcomp(s)
    n = len(s)
    for f in (0, 1, 2):
        for orf in _scan(rc[f:], f, -1):
            # reflect back to forward 1-based coords
            start = n - orf.end + 1
            end = n - orf.start + 1
            outs.append(ORF(start, end, orf.frame, -1, orf.aa_len, orf.prod))
    return sorted(outs, key=lambda o: (o.start, -o.aa_len))


def find_motifs(
    seq: str, pattern: str, *, case_insensitive: bool = True
) -> list[tuple[int, int]]:
    """Regex motif finder; returns 1-based (start,end) for each match."""
    import re

    flags = re.IGNORECASE if case_insensitive else 0
    out: list[tuple[int, int]] = []
    for m in re.finditer(pattern, seq or "", flags=flags):
        out.append((m.start() + 1, m.end()))
    return out


def find_restriction_sites(
    seq: str, enzyme_names: Iterable[str] | None = None
) -> dict[str, list[int]]:
    """
    Find cut positions (1-based) for selected enzymes using Bio.Restriction.
    Returns {enzyme: [positions...]}.
    """
    try:
        from Bio.Restriction import AllEnzymes, RestrictionBatch
        from Bio.Seq import Seq
    except Exception:
        return {}

    seq_obj = Seq(seq or "")
    enzymes = (
        RestrictionBatch(list(AllEnzymes))
        if not enzyme_names
        else RestrictionBatch(enzyme_names)
    )
    analysis = enzymes.search(seq_obj)
    # analysis is {Enz: [positions]}
    return {str(k): v for k, v in analysis.items() if v}
