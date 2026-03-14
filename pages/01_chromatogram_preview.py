from __future__ import annotations

from dataclasses import replace
from typing import cast

import streamlit as st

from abi_sauce.chromatogram import build_chromatogram_view
from abi_sauce.chromatogram_figure import build_chromatogram_figure
from abi_sauce.export import to_fasta
from abi_sauce.models import SequenceOrientation
from abi_sauce.orientation import (
    orient_left_right_values,
    orient_trim_config_for_display,
    raw_trim_config_from_display,
)
from abi_sauce.services.batch import apply_trim_configs
from abi_sauce.trim_state import (
    BatchTrimState,
    apply_submitted_trim_config,
    build_record_annotations,
    resolve_active_trim_config,
    resolve_batch_trim_inputs,
)
from abi_sauce.trimming import TrimConfig
from abi_sauce.upload_state import get_active_parsed_batch, update_active_parsed_record
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
_SAMPLE_SELECTED_RECORD_WIDGET_KEY = "sample_viewer.selected_record_widget"
_SAMPLE_TRIM_FORM_SCOPE_KEY = "sample_viewer.trim_form_scope"
_SAMPLE_TRIM_FORM_RECORD_NAME_KEY = "sample_viewer.trim_form_record_name"
_SAMPLE_TRIM_FORM_CONFIG_KEY = "sample_viewer.trim_form_config"
_SAMPLE_TRIM_FORM_ORIENTATION_KEY = "sample_viewer.trim_form_orientation"
_SAMPLE_ORIENTATION_WIDGET_KEY = "sample_viewer.orientation_is_reverse_complement"
_SAMPLE_ORIENTATION_RECORD_NAME_KEY = "sample_viewer.orientation_record_name"

st.set_page_config(page_title="Sample Viewer", layout="wide")


def _load_trim_config_into_form(
    config: TrimConfig,
    *,
    orientation: SequenceOrientation,
) -> None:
    display_trim_config = orient_trim_config_for_display(config, orientation)
    st.session_state[_SAMPLE_TRIM_LEFT_KEY] = display_trim_config.left_trim
    st.session_state[_SAMPLE_TRIM_RIGHT_KEY] = display_trim_config.right_trim
    st.session_state[_SAMPLE_TRIM_MIN_LENGTH_KEY] = config.min_length
    st.session_state[_SAMPLE_TRIM_QUALITY_ENABLED_KEY] = config.quality_trim_enabled
    st.session_state[_SAMPLE_TRIM_ERROR_CUTOFF_KEY] = config.error_probability_cutoff


def _trim_config_from_form(
    *,
    orientation: SequenceOrientation,
) -> TrimConfig:
    display_trim_config = TrimConfig(
        left_trim=int(st.session_state[_SAMPLE_TRIM_LEFT_KEY]),
        right_trim=int(st.session_state[_SAMPLE_TRIM_RIGHT_KEY]),
        min_length=int(st.session_state[_SAMPLE_TRIM_MIN_LENGTH_KEY]),
        quality_trim_enabled=bool(st.session_state[_SAMPLE_TRIM_QUALITY_ENABLED_KEY]),
        error_probability_cutoff=float(st.session_state[_SAMPLE_TRIM_ERROR_CUTOFF_KEY]),
    )
    return raw_trim_config_from_display(display_trim_config, orientation)


def _resolve_sample_form_trim_config(
    trim_state: BatchTrimState,
    *,
    selected_record_name: str,
) -> TrimConfig:
    record_trim_config = trim_state.trim_configs_by_record.get(selected_record_name)
    if record_trim_config is not None:
        return record_trim_config

    if trim_state.trim_scope == "all" and trim_state.global_trim_config is not None:
        return trim_state.global_trim_config

    return TrimConfig()


def _apply_trim_form() -> None:
    parsed_batch = get_active_parsed_batch(st.session_state)
    selected_record_name = cast(str, get_selected_record_name(st.session_state))
    if parsed_batch is None:
        return

    record = parsed_batch.parsed_records.get(selected_record_name)
    if record is None:
        return

    current_trim_state = get_batch_trim_state(st.session_state)
    updated_trim_state = apply_submitted_trim_config(
        BatchTrimState(
            trim_scope="selected",
            global_trim_config=current_trim_state.global_trim_config,
            trim_configs_by_record=dict(current_trim_state.trim_configs_by_record),
        ),
        selected_record_name=selected_record_name,
        submitted_trim_config=_trim_config_from_form(orientation=record.orientation),
    )
    set_batch_trim_state(st.session_state, updated_trim_state)


def _apply_orientation_toggle() -> None:
    parsed_batch = get_active_parsed_batch(st.session_state)
    selected_record_name = get_selected_record_name(st.session_state)
    if parsed_batch is None or selected_record_name is None:
        return

    record = parsed_batch.parsed_records.get(selected_record_name)
    if record is None:
        return

    next_orientation: SequenceOrientation = (
        "reverse_complement"
        if bool(st.session_state.get(_SAMPLE_ORIENTATION_WIDGET_KEY))
        else "forward"
    )
    if record.orientation == next_orientation:
        return

    update_active_parsed_record(
        st.session_state,
        source_filename=selected_record_name,
        record=replace(record, orientation=next_orientation),
    )


parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info(
        "Load one or more .ab1 files from the workspace sidebar to "
        "inspect individual samples."
    )
    st.stop()

parsed_records = parsed_batch.parsed_records
parse_errors = parsed_batch.parse_errors


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
sample_col, orientation_col = st.columns([4, 2])
with sample_col:
    selected_record_name = cast(
        str,
        st.selectbox(
            "Sample",
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

orientation_is_reverse_complement = record.orientation == "reverse_complement"
should_refresh_orientation_widget = (
    st.session_state.get(_SAMPLE_ORIENTATION_RECORD_NAME_KEY) != selected_record_name
    or st.session_state.get(_SAMPLE_ORIENTATION_WIDGET_KEY)
    != orientation_is_reverse_complement
)
if should_refresh_orientation_widget:
    st.session_state[_SAMPLE_ORIENTATION_WIDGET_KEY] = orientation_is_reverse_complement
    st.session_state[_SAMPLE_ORIENTATION_RECORD_NAME_KEY] = selected_record_name

record_annotations = build_record_annotations(
    parsed_records.keys(),
    trim_state.trim_configs_by_record,
)
active_form_trim_config = _resolve_sample_form_trim_config(
    trim_state,
    selected_record_name=selected_record_name,
)
should_refresh_trim_form = (
    st.session_state.get(_SAMPLE_TRIM_FORM_SCOPE_KEY) != trim_state.trim_scope
    or st.session_state.get(_SAMPLE_TRIM_FORM_RECORD_NAME_KEY) != selected_record_name
    or st.session_state.get(_SAMPLE_TRIM_FORM_CONFIG_KEY) != active_form_trim_config
    or st.session_state.get(_SAMPLE_TRIM_FORM_ORIENTATION_KEY) != record.orientation
)
if should_refresh_trim_form:
    _load_trim_config_into_form(
        active_form_trim_config,
        orientation=record.orientation,
    )
    st.session_state[_SAMPLE_TRIM_FORM_SCOPE_KEY] = trim_state.trim_scope
    st.session_state[_SAMPLE_TRIM_FORM_RECORD_NAME_KEY] = selected_record_name
    st.session_state[_SAMPLE_TRIM_FORM_CONFIG_KEY] = active_form_trim_config
    st.session_state[_SAMPLE_TRIM_FORM_ORIENTATION_KEY] = record.orientation
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
display_bases_removed_left, display_bases_removed_right = orient_left_right_values(
    trim_result.bases_removed_left,
    trim_result.bases_removed_right,
    record.orientation,
)
chromatogram_view = build_chromatogram_view(record, trim_result)

st.toggle(
    "Reverse-complement",
    key=_SAMPLE_ORIENTATION_WIDGET_KEY,
    help=(
        "Persist this sample in reverse-complement orientation for "
        "sequence preview, export output, and chromatogram display."
    ),
    on_change=_apply_orientation_toggle,
)

if not chromatogram_view.is_renderable:
    st.warning("Selected record does not have enough trace data to render a chart.")
    st.write(
        {
            "render_failure_reason": chromatogram_view.render_failure_reason,
        }
    )
else:
    figure = build_chromatogram_figure(chromatogram_view)
    figure.update_layout(
        height=500,
        margin={"l": 24, "r": 24, "t": 24, "b": 24},
    )
    st.plotly_chart(
        figure,
        width="stretch",
        config={"scrollZoom": False},
    )

metric_col_1, metric_col_2, metric_col_3, metric_col_4, change_trim_col = st.columns(
    [3, 3, 3, 3, 2]
)
with metric_col_1:
    st.metric("Original length", trim_result.original_length)
with metric_col_2:
    st.metric("Trimmed length", trim_result.trimmed_length)
with metric_col_3:
    st.metric("Bases removed left", display_bases_removed_left)
with metric_col_4:
    st.metric("Bases removed right", display_bases_removed_right)
with change_trim_col:
    with st.popover(
        "Modify trimming",
        key="modify_trim",
        width="stretch",
        icon=":material/discover_tune:",
    ):
        if record.orientation == "reverse_complement":
            st.caption(
                "Trim inputs are interpreted relative to the displayed reverse-complement view."
            )

        trim_left_col, trim_right_col, trim_min_length_col = st.columns(3)
        with trim_left_col:
            st.number_input(
                "Left trim",
                min_value=0,
                step=1,
                key=_SAMPLE_TRIM_LEFT_KEY,
                on_change=_apply_trim_form,
            )
        with trim_right_col:
            st.number_input(
                "Right trim",
                min_value=0,
                step=1,
                key=_SAMPLE_TRIM_RIGHT_KEY,
                on_change=_apply_trim_form,
            )
        with trim_min_length_col:
            st.number_input(
                "Minimum length",
                min_value=0,
                step=1,
                key=_SAMPLE_TRIM_MIN_LENGTH_KEY,
                on_change=_apply_trim_form,
            )

        quality_trim_enabled = st.checkbox(
            "Enable Mott quality trimming",
            key=_SAMPLE_TRIM_QUALITY_ENABLED_KEY,
            on_change=_apply_trim_form,
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
            on_change=_apply_trim_form,
        )

trimmed_fasta = to_fasta(trim_result.record, line_width=None)
st.code(trimmed_fasta, wrap_lines=True)

debug = False
if debug:
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

    with st.expander("Chromatogram debug info"):
        st.write(
            {
                "trace_length": chromatogram_view.trace_length,
                "channels": [
                    f"{channel.data_key}:{channel.base}"
                    for channel in chromatogram_view.channels
                ],
                "base_calls": len(chromatogram_view.base_calls),
                "quality_segments": len(chromatogram_view.quality_segments),
                "left_trim_boundary": chromatogram_view.trim_boundaries.left,
                "right_trim_boundary": chromatogram_view.trim_boundaries.right,
            }
        )

    with st.expander("Batch/session info"):
        st.write(
            {
                "uploaded_files": len(parsed_batch.uploads),
                "parsed_files": len(parsed_records),
                "failed_files": len(parse_errors),
            }
        )
        if parse_errors:
            for filename, message in parse_errors.items():
                st.error(f"{filename}: {message}")

    raw_col, trimmed_col = st.columns(2)
    with raw_col:
        st.subheader("Raw sequence")
        st.code(record.sequence[:1000] or "<empty>")
    with trimmed_col:
        st.subheader("Trimmed sequence")
        st.code(trim_result.record.sequence[:1000] or "<empty>")
