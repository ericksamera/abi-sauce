from __future__ import annotations

from abi_sauce.reference_alignment_types import AlignmentResult


def format_reference_alignment_fasta(
    result: AlignmentResult,
    *,
    reference_name: str | None = None,
    query_name: str | None = None,
) -> str:
    """Render one reference-guided alignment as two gapped FASTA entries."""
    resolved_reference_name = (
        reference_name.strip()
        if isinstance(reference_name, str) and reference_name.strip()
        else result.reference_name
    )
    resolved_query_name = (
        query_name.strip()
        if isinstance(query_name, str) and query_name.strip()
        else f"{result.sample_name}__{result.strand}"
    )

    return (
        f">{resolved_reference_name}\n"
        f"{result.aligned_reference}\n"
        f">{resolved_query_name}\n"
        f"{result.aligned_query}\n"
    )
