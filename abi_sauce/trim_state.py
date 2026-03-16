from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Literal

from abi_sauce.trimming import TrimConfig

TrimScope = Literal["all", "selected"]

DEFAULT_BATCH_TRIM_CONFIG = TrimConfig(
    quality_trim_enabled=True,
    error_probability_cutoff=0.01,
    min_length=25,
)


@dataclass(frozen=True, slots=True)
class BatchTrimState:
    """Pure UI-level trimming state for batch pages."""

    trim_scope: TrimScope = "all"
    global_trim_config: TrimConfig | None = DEFAULT_BATCH_TRIM_CONFIG
    trim_configs_by_record: Mapping[str, TrimConfig] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ResolvedBatchTrimInputs:
    """Effective trim inputs to pass into the batch service."""

    default_trim_config: TrimConfig | None = None
    trim_configs_by_name: dict[str, TrimConfig] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class RecordTrimAnnotations:
    """Derived record labels and flags for the trimming UI."""

    display_labels_by_record: dict[str, str] = field(default_factory=dict)
    overridden_record_names: frozenset[str] = field(default_factory=frozenset)
    overridden_count: int = 0
    custom_trim_flags_by_record: dict[str, bool] = field(default_factory=dict)


def resolve_batch_trim_inputs(
    trim_state: BatchTrimState,
) -> ResolvedBatchTrimInputs:
    """Resolve the effective batch trim inputs for the current UI state."""
    resolved_default_trim_config = _resolved_global_trim_config(trim_state)

    if trim_state.trim_scope == "all":
        return ResolvedBatchTrimInputs(
            default_trim_config=resolved_default_trim_config,
        )

    return ResolvedBatchTrimInputs(
        default_trim_config=resolved_default_trim_config,
        trim_configs_by_name=dict(trim_state.trim_configs_by_record),
    )


def apply_submitted_trim_config(
    trim_state: BatchTrimState,
    *,
    selected_record_name: str,
    submitted_trim_config: TrimConfig,
) -> BatchTrimState:
    """Return updated trim state for one submitted form config."""
    normalized_trim_config = _normalize_trim_config(submitted_trim_config)

    if trim_state.trim_scope == "all":
        return BatchTrimState(
            trim_scope=trim_state.trim_scope,
            global_trim_config=normalized_trim_config,
            trim_configs_by_record=dict(trim_state.trim_configs_by_record),
        )

    updated_trim_configs_by_record = dict(trim_state.trim_configs_by_record)
    if normalized_trim_config is None:
        updated_trim_configs_by_record.pop(selected_record_name, None)
    else:
        updated_trim_configs_by_record[selected_record_name] = normalized_trim_config

    return BatchTrimState(
        trim_scope=trim_state.trim_scope,
        global_trim_config=trim_state.global_trim_config,
        trim_configs_by_record=updated_trim_configs_by_record,
    )


def build_record_annotations(
    record_names: Iterable[str],
    trim_configs_by_record: Mapping[str, TrimConfig],
) -> RecordTrimAnnotations:
    """Build derived labels and flags for record-level trim UI annotations."""
    record_names = tuple(record_names)
    current_record_names = frozenset(record_names)
    overridden_record_names = frozenset(
        record_name
        for record_name in trim_configs_by_record
        if record_name in current_record_names
    )

    display_labels_by_record = {
        record_name: (
            f"{record_name} *"
            if record_name in overridden_record_names
            else record_name
        )
        for record_name in record_names
    }
    custom_trim_flags_by_record = {
        record_name: record_name in overridden_record_names
        for record_name in record_names
    }

    return RecordTrimAnnotations(
        display_labels_by_record=display_labels_by_record,
        overridden_record_names=overridden_record_names,
        overridden_count=len(overridden_record_names),
        custom_trim_flags_by_record=custom_trim_flags_by_record,
    )


def resolve_active_trim_config(
    trim_state: BatchTrimState,
    *,
    selected_record_name: str,
) -> TrimConfig:
    """Resolve the trim config that should currently hydrate the form."""
    if trim_state.trim_scope == "all":
        return _resolved_global_trim_config(trim_state)

    record_trim_config = trim_state.trim_configs_by_record.get(selected_record_name)
    if record_trim_config is not None:
        return record_trim_config
    return _resolved_global_trim_config(trim_state)


def _resolved_global_trim_config(trim_state: BatchTrimState) -> TrimConfig:
    """Return the stored batch default, falling back to the built-in default."""
    return (
        trim_state.global_trim_config
        if trim_state.global_trim_config is not None
        else DEFAULT_BATCH_TRIM_CONFIG
    )


def _normalize_trim_config(config: TrimConfig) -> TrimConfig | None:
    """Normalize a no-op trim config to None for state storage."""
    return None if config == TrimConfig() else config
