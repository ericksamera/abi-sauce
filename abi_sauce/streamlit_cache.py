from __future__ import annotations

from collections.abc import Iterable
import hashlib
from typing import TypeAlias

import streamlit as st

from abi_sauce.assembly_types import AssemblyStrand, MultiAssemblyResult
from abi_sauce.assembly_state import AssemblyDefinition
from abi_sauce.assembly_trace import (
    AssemblyTraceRowSource,
    AssemblyTraceView,
    build_assembly_trace_row_source,
    build_multi_assembly_trace_view,
    build_pairwise_assembly_trace_view,
)
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.services.assembly_compute import (
    ComputedAssembly,
    compute_saved_assemblies,
)
from abi_sauce.services.batch_parse import (
    BatchSignature,
    ParsedBatch,
)
from abi_sauce.services.batch_trim import (
    PreparedBatch,
    build_prepared_batch,
    resolve_effective_trim_configs,
)
from abi_sauce.trim_state import (
    BatchTrimState,
    ResolvedBatchTrimInputs,
    resolve_batch_trim_inputs,
)
from abi_sauce.trimming import TrimConfig, TrimResult, trim_sequence_record

ParsedBatchCacheKey: TypeAlias = tuple[
    BatchSignature,
    tuple[tuple[str, str], ...],
    tuple[tuple[str, str], ...],
]
TrimInputsCacheKey: TypeAlias = tuple[
    TrimConfig | None,
    tuple[tuple[str, TrimConfig], ...],
]
PreparedBatchCacheKey: TypeAlias = tuple[
    ParsedBatchCacheKey,
    tuple[tuple[str, str], ...],
]
AssemblyDefinitionsCacheKey: TypeAlias = tuple[AssemblyDefinition, ...]
ComputedAssemblyCacheKey: TypeAlias = tuple[
    AssemblyDefinition,
    str,
    str | None,
    str | None,
    str | None,
]
AssemblyTraceRowSourceCacheKey: TypeAlias = tuple[str, str, AssemblyStrand]
TrimmedRecordCacheKey: TypeAlias = tuple[str, TrimConfig]

_PREPARED_BATCH_CACHE_VERSION = 3
_COMPUTED_ASSEMBLIES_CACHE_VERSION = 2
_ASSEMBLY_TRACE_ROW_SOURCE_CACHE_VERSION = 1
_TRIM_SEQUENCE_RECORD_CACHE_VERSION = 1
_ASSEMBLY_TRACE_VIEW_CACHE_VERSION = 2


def build_parsed_batch_cache_key(
    parsed_batch: ParsedBatch,
) -> ParsedBatchCacheKey:
    """Return a stable cache key for one parsed batch."""
    return (
        parsed_batch.signature,
        tuple(
            sorted(
                (
                    source_filename,
                    _sequence_record_cache_digest(record),
                )
                for source_filename, record in parsed_batch.parsed_records.items()
            )
        ),
        tuple(sorted(parsed_batch.parse_errors.items())),
    )


def build_trim_inputs_cache_key(
    resolved_trim_inputs: ResolvedBatchTrimInputs,
) -> TrimInputsCacheKey:
    """Return a stable hashable cache key for resolved batch trim inputs."""
    return (
        resolved_trim_inputs.default_trim_config,
        tuple(sorted(resolved_trim_inputs.trim_configs_by_name.items())),
    )


def build_prepared_batch_cache_key(
    prepared_batch: PreparedBatch,
) -> PreparedBatchCacheKey:
    """Return a stable cache key for one prepared batch."""
    return (
        build_parsed_batch_cache_key(
            ParsedBatch(
                uploads=prepared_batch.uploads,
                parsed_records=prepared_batch.parsed_records,
                parse_errors=prepared_batch.parse_errors,
                signature=prepared_batch.signature,
            )
        ),
        tuple(
            sorted(
                (
                    source_filename,
                    _trim_result_cache_digest(trim_result),
                )
                for source_filename, trim_result in prepared_batch.trim_results.items()
            )
        ),
    )


def build_assembly_definitions_cache_key(
    definitions: Iterable[AssemblyDefinition],
) -> AssemblyDefinitionsCacheKey:
    """Return a stable cache key for saved assembly definitions."""
    return tuple(definitions)


def build_computed_assembly_cache_key(
    computed_assembly: ComputedAssembly,
) -> ComputedAssemblyCacheKey:
    """Return a stable cache key for one computed assembly result."""
    return (
        computed_assembly.definition,
        computed_assembly.status,
        computed_assembly.status_reason,
        (
            None
            if computed_assembly.result is None
            else _stable_repr_digest(computed_assembly.result)
        ),
        (
            None
            if computed_assembly.consensus_record is None
            else _sequence_record_cache_digest(computed_assembly.consensus_record)
        ),
    )


def build_assembly_trace_row_source_cache_key(
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    *,
    strand: AssemblyStrand,
) -> AssemblyTraceRowSourceCacheKey:
    """Return a stable cache key for one aligned-trace row source."""
    return (
        _sequence_record_cache_digest(raw_record),
        _trim_result_cache_digest(trim_result),
        strand,
    )


def build_trimmed_record_cache_key(
    parsed_record: SequenceRecord,
    trim_config: TrimConfig,
) -> TrimmedRecordCacheKey:
    """Return a stable cache key for one record/config trim computation."""
    return (
        _sequence_record_cache_digest(parsed_record),
        trim_config,
    )


@st.cache_data(show_spinner=False, max_entries=4096)
def _trim_sequence_record_cached(
    *,
    cache_version: int,
    trimmed_record_key: TrimmedRecordCacheKey,
    _parsed_record: SequenceRecord,
    _trim_config: TrimConfig,
) -> TrimResult:
    """Return one cached trim result for a single record/config pair."""
    return trim_sequence_record(_parsed_record, _trim_config)


def trim_sequence_record_for_cache(
    parsed_record: SequenceRecord,
    trim_config: TrimConfig,
) -> TrimResult:
    """Return one cached trim result for a single record/config pair."""
    return _trim_sequence_record_cached(
        cache_version=_TRIM_SEQUENCE_RECORD_CACHE_VERSION,
        trimmed_record_key=build_trimmed_record_cache_key(
            parsed_record,
            trim_config,
        ),
        _parsed_record=parsed_record,
        _trim_config=trim_config,
    )


@st.cache_data(show_spinner=False, max_entries=256)
def _build_assembly_trace_row_source_cached(
    *,
    cache_version: int,
    row_source_key: AssemblyTraceRowSourceCacheKey,
    _raw_record: SequenceRecord,
    _trim_result: TrimResult,
) -> AssemblyTraceRowSource:
    """Return a cached aligned-trace row source for one member."""
    _record_digest, _trim_digest, strand = row_source_key
    return build_assembly_trace_row_source(
        raw_record=_raw_record,
        trim_result=_trim_result,
        strand=strand,
    )


def build_assembly_trace_row_source_for_member(
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    *,
    strand: AssemblyStrand,
) -> AssemblyTraceRowSource:
    """Return a cached aligned-trace row source for one assembly member."""
    return _build_assembly_trace_row_source_cached(
        cache_version=_ASSEMBLY_TRACE_ROW_SOURCE_CACHE_VERSION,
        row_source_key=build_assembly_trace_row_source_cache_key(
            raw_record,
            trim_result,
            strand=strand,
        ),
        _raw_record=raw_record,
        _trim_result=trim_result,
    )


@st.cache_data(show_spinner=False, max_entries=16)
def _prepare_batch_cached(
    *,
    cache_version: int,
    parsed_batch_key: ParsedBatchCacheKey,
    trim_inputs_key: TrimInputsCacheKey,
    _parsed_batch: ParsedBatch,
    _resolved_trim_inputs: ResolvedBatchTrimInputs,
) -> PreparedBatch:
    """Return one prepared batch for a parsed batch plus resolved trim inputs."""
    effective_trim_configs_by_name = resolve_effective_trim_configs(
        _parsed_batch,
        default_trim_config=_resolved_trim_inputs.default_trim_config,
        trim_configs_by_name=_resolved_trim_inputs.trim_configs_by_name,
    )
    trim_results = {
        source_filename: trim_sequence_record_for_cache(
            _parsed_batch.parsed_records[source_filename],
            effective_trim_configs_by_name[source_filename],
        )
        for source_filename in _parsed_batch.parsed_records
    }
    return build_prepared_batch(
        _parsed_batch,
        trim_results=trim_results,
    )


def prepare_batch_for_trim_inputs(
    parsed_batch: ParsedBatch,
    resolved_trim_inputs: ResolvedBatchTrimInputs,
) -> PreparedBatch:
    """Return a cached prepared batch for the current parsed batch and trim inputs."""
    return _prepare_batch_cached(
        cache_version=_PREPARED_BATCH_CACHE_VERSION,
        parsed_batch_key=build_parsed_batch_cache_key(parsed_batch),
        trim_inputs_key=build_trim_inputs_cache_key(resolved_trim_inputs),
        _parsed_batch=parsed_batch,
        _resolved_trim_inputs=resolved_trim_inputs,
    )


def prepare_batch_for_trim_state(
    parsed_batch: ParsedBatch,
    trim_state: BatchTrimState,
) -> PreparedBatch:
    """Resolve batch trim state and return the cached prepared batch."""
    return prepare_batch_for_trim_inputs(
        parsed_batch,
        resolve_batch_trim_inputs(trim_state),
    )


@st.cache_data(show_spinner=False, max_entries=64)
def _compute_saved_assemblies_cached(
    *,
    cache_version: int,
    prepared_batch_key: PreparedBatchCacheKey,
    definitions_key: AssemblyDefinitionsCacheKey,
    _prepared_batch: PreparedBatch,
    _definitions: tuple[AssemblyDefinition, ...],
) -> dict[str, ComputedAssembly]:
    """Return saved assemblies recomputed against one prepared batch."""
    return compute_saved_assemblies(_prepared_batch, _definitions)


def compute_saved_assemblies_for_definitions(
    prepared_batch: PreparedBatch,
    definitions: Iterable[AssemblyDefinition],
) -> dict[str, ComputedAssembly]:
    """Return cached saved-assembly results for the current prepared batch."""
    definitions_tuple = tuple(definitions)
    return _compute_saved_assemblies_cached(
        cache_version=_COMPUTED_ASSEMBLIES_CACHE_VERSION,
        prepared_batch_key=build_prepared_batch_cache_key(prepared_batch),
        definitions_key=build_assembly_definitions_cache_key(definitions_tuple),
        _prepared_batch=prepared_batch,
        _definitions=definitions_tuple,
    )


@st.cache_data(show_spinner=False, max_entries=64)
def _build_selected_assembly_trace_view_cached(
    *,
    cache_version: int,
    prepared_batch_key: PreparedBatchCacheKey,
    computed_assembly_key: ComputedAssemblyCacheKey,
    cell_width: float,
    samples_per_cell: int | None,
    trace_row_height: float,
    _prepared_batch: PreparedBatch,
    _computed_assembly: ComputedAssembly,
) -> AssemblyTraceView | None:
    """Return a cached aligned trace view for one selected computed assembly."""
    result = _computed_assembly.result
    if result is None:
        return None

    if isinstance(result, MultiAssemblyResult):
        row_sources_by_member_index = {
            member.member_index: build_assembly_trace_row_source_for_member(
                _prepared_batch.parsed_records[member.source_filename],
                _prepared_batch.trim_results[member.source_filename],
                strand=member.chosen_orientation,
            )
            for member in result.members
            if member.included
        }
        return build_multi_assembly_trace_view(
            result=result,
            raw_records_by_source_filename=_prepared_batch.parsed_records,
            trim_results_by_source_filename=_prepared_batch.trim_results,
            cell_width=cell_width,
            samples_per_cell=samples_per_cell,
            trace_row_height=trace_row_height,
            row_sources_by_member_index=row_sources_by_member_index,
        )

    source_filenames = _computed_assembly.definition.source_filenames
    if len(source_filenames) != 2:
        return None
    left_source_filename, right_source_filename = source_filenames
    left_row_source = build_assembly_trace_row_source_for_member(
        _prepared_batch.parsed_records[left_source_filename],
        _prepared_batch.trim_results[left_source_filename],
        strand="forward",
    )
    right_row_source = build_assembly_trace_row_source_for_member(
        _prepared_batch.parsed_records[right_source_filename],
        _prepared_batch.trim_results[right_source_filename],
        strand=result.chosen_right_orientation,
    )
    return build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename=left_source_filename,
        left_raw_record=_prepared_batch.parsed_records[left_source_filename],
        left_trim_result=_prepared_batch.trim_results[left_source_filename],
        right_source_filename=right_source_filename,
        right_raw_record=_prepared_batch.parsed_records[right_source_filename],
        right_trim_result=_prepared_batch.trim_results[right_source_filename],
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
        left_row_source=left_row_source,
        right_row_source=right_row_source,
    )


def build_selected_assembly_trace_view(
    prepared_batch: PreparedBatch,
    computed_assembly: ComputedAssembly,
    *,
    cell_width: float = 1.0,
    samples_per_cell: int | None = None,
    trace_row_height: float = 3.0,
) -> AssemblyTraceView | None:
    """Return a cached aligned trace view for the selected assembly."""
    return _build_selected_assembly_trace_view_cached(
        cache_version=_ASSEMBLY_TRACE_VIEW_CACHE_VERSION,
        prepared_batch_key=build_prepared_batch_cache_key(prepared_batch),
        computed_assembly_key=build_computed_assembly_cache_key(computed_assembly),
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
        _prepared_batch=prepared_batch,
        _computed_assembly=computed_assembly,
    )


def _sequence_record_cache_digest(record: SequenceRecord) -> str:
    snapshot = (
        record.record_id,
        record.name,
        record.description,
        record.sequence,
        record.source_format,
        record.orientation,
        (
            None
            if record.qualities is None
            else tuple(int(value) for value in record.qualities)
        ),
        _trace_data_cache_snapshot(record.trace_data),
        _annotations_cache_snapshot(record.annotations),
    )
    return _stable_repr_digest(snapshot)


def _trim_result_cache_digest(trim_result: TrimResult) -> str:
    snapshot = (
        _sequence_record_cache_digest(trim_result.record),
        trim_result.original_length,
        trim_result.trimmed_length,
        trim_result.bases_removed_left,
        trim_result.bases_removed_right,
        trim_result.passed_min_length,
        trim_result.fixed_bases_removed_left,
        trim_result.fixed_bases_removed_right,
        trim_result.quality_bases_removed_left,
        trim_result.quality_bases_removed_right,
    )
    return _stable_repr_digest(snapshot)


def _trace_data_cache_snapshot(
    trace_data: TraceData | None,
) -> tuple[str | None, tuple[int, ...], tuple[tuple[str, tuple[int, ...]], ...]] | None:
    if trace_data is None:
        return None

    return (
        trace_data.channel_order,
        tuple(int(value) for value in trace_data.base_positions),
        tuple(
            (
                channel_name,
                tuple(int(value) for value in signal),
            )
            for channel_name, signal in sorted(trace_data.channels.items())
        ),
    )


def _annotations_cache_snapshot(
    annotations: dict[str, object],
) -> tuple[tuple[str, str], ...]:
    return tuple(
        sorted(
            (
                str(key),
                _stable_repr_digest(value),
            )
            for key, value in annotations.items()
        )
    )


def _stable_repr_digest(value: object) -> str:
    return hashlib.blake2b(
        repr(value).encode("utf-8"),
        digest_size=16,
    ).hexdigest()
