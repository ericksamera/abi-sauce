from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.chromatogram import build_chromatogram_view
from abi_sauce.chromatogram_figure import build_chromatogram_figure
from abi_sauce.services.batch import apply_trim_configs
from abi_sauce.trim_state import (
    BatchTrimState,
    TrimScope,
    apply_submitted_trim_config,
    build_record_annotations,
    resolve_active_trim_config,
    resolve_batch_trim_inputs,
)
from abi_sauce.trimming import TrimConfig
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import (
    get_batch_trim_state,
    get_selected_record_name,
    set_batch_trim_state,
    set_selected_record_name,
    sync_viewer_session_state,
)

_SAMPLE_TRIM_LEFT_KEY = "sample_viewer.trim_left"
_SAMPLE_TRIM_RIGHT_KEY = "sample_viewer.trim_right"
_SAMPLE_TRIM_MIN_LENGTH_KEY = "sample_viewer.trim_min_length"
_SAMPLE_TRIM_QUALITY_ENABLED_KEY = "sample_viewer.trim_quality_enabled"
_SAMPLE_TRIM_ERROR_CUTOFF_KEY = "sample_viewer.trim_error_probability_cutoff"
_SAMPLE_TRIM_SCOPE_WIDGET_KEY = "sample_viewer.trim_scope_widget"
_SAMPLE_SELECTED_RECORD_WIDGET_KEY = "sample_viewer.selected_record_widget"
_SAMPLE_TRIM_FORM_SCOPE_KEY = "sample_viewer.trim_form_scope"
_SAMPLE_TRIM_FORM_RECORD_NAME_KEY = "sample_viewer.trim_form_record_name"

st.set_page_config(page_title="Sample Viewer", layout="wide")
st.title("Sample Viewer")
st.caption(
    "Inspect one selected sample from the shared active batch, edit shared trim state, and visualize the chromatogram."
)


def _load_trim_config_into_form(config: TrimConfig) -> None:
    st.session_state[_SAMPLE_TRIM_LEFT_KEY] = config.left_trim
    st.session_state[_SAMPLE_TRIM_RIGHT_KEY] = config.right_trim
    st.session_state[_SAMPLE_TRIM_MIN_LENGTH_KEY] = config.min_length
    st.session_state[_SAMPLE_TRIM_QUALITY_ENABLED_KEY] = config.quality_trim_enabled
    st.session_state[_SAMPLE_TRIM_ERROR_CUTOFF_KEY] = config.error_probability_cutoff


def _trim_config_from_form() -> TrimConfig:
    return TrimConfig(
        left_trim=int(st.session_state[_SAMPLE_TRIM_LEFT_KEY]),
        right_trim=int(st.session_state[_SAMPLE_TRIM_RIGHT_KEY]),
        min_length=int(st.session_state[_SAMPLE_TRIM_MIN_LENGTH_KEY]),
        quality_trim_enabled=bool(st.session_state[_SAMPLE_TRIM_QUALITY_ENABLED_KEY]),
        error_probability_cutoff=float(st.session_state[_SAMPLE_TRIM_ERROR_CUTOFF_KEY]),
    )


def _apply_trim_form() -> None:
    updated_trim_state = apply_submitted_trim_config(
        get_batch_trim_state(st.session_state),
        selected_record_name=cast(str, get_selected_record_name(st.session_state)),
        submitted_trim_config=_trim_config_from_form(),
    )
    set_batch_trim_state(st.session_state, updated_trim_state)


parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info(
        "Load one or more .ab1 files from the workspace sidebar to inspect individual samples."
    )
    st.stop()

parsed_records = parsed_batch.parsed_records
parse_errors = parsed_batch.parse_errors

st.write(
    {
        "uploaded_files": len(parsed_batch.uploads),
        "parsed_files": len(parsed_records),
        "failed_files": len(parse_errors),
    }
)

if parse_errors:
    with st.expander("Files that failed to parse"):
        for filename, message in parse_errors.items():
            st.error(f"{filename}: {message}")

if not parsed_records:
    st.error("No ABI files could be parsed from this batch.")
    st.stop()

viewer_state = sync_viewer_session_state(
    st.session_state,
    batch_signature=parsed_batch.signature,
    parsed_record_names=tuple(parsed_records),
)
trim_state = viewer_state.trim_state
record_names = list(parsed_records.keys())
selected_record_name = viewer_state.selected_record_name or record_names[0]

record_annotations = build_record_annotations(
    parsed_records.keys(),
    trim_state.trim_configs_by_record,
)
selected_record_name = cast(
    str,
    st.selectbox(
        "Selected record",
        options=record_names,
        index=record_names.index(selected_record_name),
        format_func=lambda filename: record_annotations.display_labels_by_record[
            filename
        ],
        key=_SAMPLE_SELECTED_RECORD_WIDGET_KEY,
    ),
)
set_selected_record_name(st.session_state, selected_record_name)
record = parsed_records[selected_record_name]

trim_scope = cast(
    TrimScope,
    st.radio(
        "Shared trim application mode",
        options=["selected", "all"],
        format_func=lambda value: (
            "Selected record override" if value == "selected" else "All parsed records"
        ),
        index=0 if trim_state.trim_scope == "selected" else 1,
        horizontal=True,
        key=_SAMPLE_TRIM_SCOPE_WIDGET_KEY,
    ),
)
if trim_scope != trim_state.trim_scope:
    trim_state = BatchTrimState(
        trim_scope=trim_scope,
        global_trim_config=trim_state.global_trim_config,
        trim_configs_by_record=dict(trim_state.trim_configs_by_record),
    )
    set_batch_trim_state(st.session_state, trim_state)

record_annotations = build_record_annotations(
    parsed_records.keys(),
    trim_state.trim_configs_by_record,
)
active_form_trim_config = resolve_active_trim_config(
    trim_state,
    selected_record_name=selected_record_name,
)
should_refresh_trim_form = st.session_state.get(
    _SAMPLE_TRIM_FORM_SCOPE_KEY
) != trim_scope or (
    trim_scope == "selected"
    and st.session_state.get(_SAMPLE_TRIM_FORM_RECORD_NAME_KEY) != selected_record_name
)
if should_refresh_trim_form:
    _load_trim_config_into_form(active_form_trim_config)
    st.session_state[_SAMPLE_TRIM_FORM_SCOPE_KEY] = trim_scope
    st.session_state[_SAMPLE_TRIM_FORM_RECORD_NAME_KEY] = selected_record_name

st.caption(
    "In all-record mode, this form edits the shared batch trim config. In selected-record mode, applying a no-op config removes the override."
)
st.caption(
    "* = custom per-sequence trim override "
    f"({record_annotations.overridden_count}/{len(parsed_records)} records marked)"
)
with st.form("sample_viewer_trim_form"):
    trim_left_col, trim_right_col, trim_min_length_col = st.columns(3)
    with trim_left_col:
        st.number_input(
            "Left trim",
            min_value=0,
            step=1,
            key=_SAMPLE_TRIM_LEFT_KEY,
        )
    with trim_right_col:
        st.number_input(
            "Right trim",
            min_value=0,
            step=1,
            key=_SAMPLE_TRIM_RIGHT_KEY,
        )
    with trim_min_length_col:
        st.number_input(
            "Minimum length",
            min_value=0,
            step=1,
            key=_SAMPLE_TRIM_MIN_LENGTH_KEY,
        )

    quality_trim_enabled = st.checkbox(
        "Enable Mott quality trimming",
        key=_SAMPLE_TRIM_QUALITY_ENABLED_KEY,
    )
    st.number_input(
        "Mott max acceptable error probability",
        min_value=0.0,
        max_value=1.0,
        step=0.0001,
        format="%.4f",
        key=_SAMPLE_TRIM_ERROR_CUTOFF_KEY,
        disabled=not quality_trim_enabled,
        help=(
            "Lower values trim more aggressively. Examples: "
            "Q20 = 0.01, Q30 = 0.001, Q40 = 0.0001."
        ),
    )

    st.form_submit_button("Apply shared trim", on_click=_apply_trim_form)

trim_state = get_batch_trim_state(st.session_state)
resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)
prepared_batch = apply_trim_configs(
    parsed_batch,
    default_trim_config=resolved_trim_inputs.default_trim_config,
    trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
)
trim_result = prepared_batch.trim_results[selected_record_name]
effective_trim_config = resolve_active_trim_config(
    trim_state,
    selected_record_name=selected_record_name,
)
chromatogram_view = build_chromatogram_view(record, trim_result)

st.subheader("Selected record detail")
metric_col_1, metric_col_2, metric_col_3, metric_col_4 = st.columns(4)
with metric_col_1:
    st.metric("Original length", trim_result.original_length)
with metric_col_2:
    st.metric("Trimmed length", trim_result.trimmed_length)
with metric_col_3:
    st.metric("Bases removed left", trim_result.bases_removed_left)
with metric_col_4:
    st.metric("Bases removed right", trim_result.bases_removed_right)

st.write(
    {
        "selected_record": selected_record_name,
        "record_id": record.record_id,
        "name": record.name,
        "has_trace_data": record.trace_data is not None,
        "quality_count": 0 if record.qualities is None else len(record.qualities),
        "trim_scope": trim_state.trim_scope,
        "quality_trim_enabled": effective_trim_config.quality_trim_enabled,
        "error_probability_cutoff": effective_trim_config.error_probability_cutoff,
        "custom_trim_override": record_annotations.custom_trim_flags_by_record.get(
            selected_record_name,
            False,
        ),
        "passed_min_length": trim_result.passed_min_length,
    }
)

if not trim_result.passed_min_length:
    st.warning("Trimmed sequence did not meet the minimum length.")

st.subheader("Chromatogram")
if not chromatogram_view.is_renderable:
    st.warning("Selected record does not have enough trace data to render a chart.")
    st.write(
        {
            "render_failure_reason": chromatogram_view.render_failure_reason,
        }
    )
else:
    st.plotly_chart(
        build_chromatogram_figure(chromatogram_view),
        width="stretch",
    )
    with st.expander("Chromatogram debug info"):
        st.write(
            {
                "trace_length": chromatogram_view.trace_length,
                "channels": [
                    f"{channel.data_key}:{channel.base}"
                    for channel in chromatogram_view.channels
                ],
                "base_calls": len(chromatogram_view.base_calls),
                "quality_points": len(chromatogram_view.quality_points),
                "left_trim_boundary": chromatogram_view.trim_boundaries.left,
                "right_trim_boundary": chromatogram_view.trim_boundaries.right,
            }
        )

raw_col, trimmed_col = st.columns(2)
with raw_col:
    st.subheader("Raw sequence")
    st.code(record.sequence[:1000] or "<empty>")
with trimmed_col:
    st.subheader("Trimmed sequence")
    st.code(trim_result.record.sequence[:1000] or "<empty>")
