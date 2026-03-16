from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from abi_sauce.batch import (
    BatchExportPolicy,
    BatchSummary,
    build_batch_export_policy,
    build_batch_summary,
)
from abi_sauce.models import SequenceRecord, SequenceUpload
from abi_sauce.services.batch_parse import BatchSignature, ParsedBatch
from abi_sauce.trimming import TrimConfig, TrimResult, trim_sequence_record


@dataclass(frozen=True, slots=True)
class PreparedBatch:
    """Parsed batch plus derived trim/export state."""

    uploads: tuple[SequenceUpload, ...]
    parsed_records: dict[str, SequenceRecord]
    parse_errors: dict[str, str]
    signature: BatchSignature
    trim_results: dict[str, TrimResult]
    batch_summary: BatchSummary
    batch_export_policy: BatchExportPolicy


def apply_trim_config(
    parsed_batch: ParsedBatch,
    trim_config: TrimConfig,
) -> PreparedBatch:
    """Apply one trim config across an entire parsed batch."""
    return apply_trim_configs(
        parsed_batch,
        default_trim_config=trim_config,
    )


def apply_trim_configs(
    parsed_batch: ParsedBatch,
    *,
    default_trim_config: TrimConfig | None = None,
    trim_configs_by_name: Mapping[str, TrimConfig] | None = None,
) -> PreparedBatch:
    """Apply trimming across a parsed batch with optional per-record configs."""
    effective_trim_configs_by_name = resolve_effective_trim_configs(
        parsed_batch,
        default_trim_config=default_trim_config,
        trim_configs_by_name=trim_configs_by_name,
    )
    trim_results = {
        name: trim_sequence_record(
            parsed_batch.parsed_records[name],
            effective_trim_configs_by_name[name],
        )
        for name in parsed_batch.parsed_records
    }
    return build_prepared_batch(
        parsed_batch,
        trim_results=trim_results,
    )


def resolve_effective_trim_configs(
    parsed_batch: ParsedBatch,
    *,
    default_trim_config: TrimConfig | None = None,
    trim_configs_by_name: Mapping[str, TrimConfig] | None = None,
) -> dict[str, TrimConfig]:
    """Resolve the effective trim config for each parsed record in the batch."""
    resolved_default_trim_config = (
        TrimConfig() if default_trim_config is None else default_trim_config
    )
    resolved_trim_configs_by_name = dict(trim_configs_by_name or {})

    return {
        name: resolved_trim_configs_by_name.get(name, resolved_default_trim_config)
        for name in parsed_batch.parsed_records
    }


def build_prepared_batch(
    parsed_batch: ParsedBatch,
    *,
    trim_results: Mapping[str, TrimResult],
) -> PreparedBatch:
    """Assemble one PreparedBatch from parsed records plus trimmed results."""
    resolved_trim_results = dict(trim_results)

    return PreparedBatch(
        uploads=parsed_batch.uploads,
        parsed_records=parsed_batch.parsed_records,
        parse_errors=parsed_batch.parse_errors,
        signature=parsed_batch.signature,
        trim_results=resolved_trim_results,
        batch_summary=build_batch_summary(
            uploads=parsed_batch.uploads,
            trim_results=resolved_trim_results,
            parse_errors=parsed_batch.parse_errors,
        ),
        batch_export_policy=build_batch_export_policy(
            trim_results=resolved_trim_results
        ),
    )


__all__ = [
    "PreparedBatch",
    "apply_trim_config",
    "apply_trim_configs",
    "build_prepared_batch",
    "resolve_effective_trim_configs",
]
