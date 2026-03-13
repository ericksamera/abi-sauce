from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.batch import ExportFormat
from abi_sauce.exceptions import ExportError
from abi_sauce.services.batch import apply_trim_configs, prepare_batch_download
from abi_sauce.trim_state import (
    BatchTrimState,
    TrimScope,
    build_record_annotations,
    resolve_batch_trim_inputs,
)
from abi_sauce.trimming import TrimConfig
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import get_batch_trim_state, sync_viewer_session_state

_BATCH_TRIM_LEFT_KEY = "batch_viewer.trim_left"
_BATCH_TRIM_RIGHT_KEY = "batch_viewer.trim_right"
_BATCH_TRIM_MIN_LENGTH_KEY = "batch_viewer.trim_min_length"
_BATCH_TRIM_QUALITY_ENABLED_KEY = "batch_viewer.trim_quality_enabled"
_BATCH_TRIM_ERROR_CUTOFF_KEY = "batch_viewer.trim_error_probability_cutoff"
_BATCH_TRIM_SCOPE_WIDGET_KEY = "batch_viewer.trim_scope_widget"
_BATCH_TRIM_FORM_SCOPE_KEY = "batch_viewer.trim_form_scope"
_BATCH_TRIM_FORM_GLOBAL_CONFIG_KEY = "batch_viewer.trim_form_global_config"

st.set_page_config(page_title="ABI Sauce", layout="wide")
st.title("Batch Viewer")
st.caption(
    "Manage shared batch-wide trim settings, review all parsed records, and prepare exports from the shared active batch."
)


def _load_trim_config_into_form(config: TrimConfig) -> None:
    st.session_state[_BATCH_TRIM_LEFT_KEY] = config.left_trim
    st.session_state[_BATCH_TRIM_RIGHT_KEY] = config.right_trim
    st.session_state[_BATCH_TRIM_MIN_LENGTH_KEY] = config.min_length
    st.session_state[_BATCH_TRIM_QUALITY_ENABLED_KEY] = config.quality_trim_enabled
    st.session_state[_BATCH_TRIM_ERROR_CUTOFF_KEY] = config.error_probability_cutoff


def _trim_config_from_form() -> TrimConfig:
    return TrimConfig(
        left_trim=int(st.session_state[_BATCH_TRIM_LEFT_KEY]),
        right_trim=int(st.session_state[_BATCH_TRIM_RIGHT_KEY]),
        min_length=int(st.session_state[_BATCH_TRIM_MIN_LENGTH_KEY]),
        quality_trim_enabled=bool(st.session_state[_BATCH_TRIM_QUALITY_ENABLED_KEY]),
        error_probability_cutoff=float(st.session_state[_BATCH_TRIM_ERROR_CUTOFF_KEY]),
    )


def _normalize_trim_config(config: TrimConfig) -> TrimConfig | None:
    return None if config == TrimConfig() else config


def _apply_global_trim_form() -> None:
    current_trim_state = get_batch_trim_state(st.session_state)
    updated_trim_state = BatchTrimState(
        trim_scope=current_trim_state.trim_scope,
        global_trim_config=_normalize_trim_config(_trim_config_from_form()),
        trim_configs_by_record=dict(current_trim_state.trim_configs_by_record),
    )
    from abi_sauce.viewer_state import set_batch_trim_state

    set_batch_trim_state(st.session_state, updated_trim_state)


parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info("Load one or more .ab1 files from the workspace sidebar to begin.")
    st.stop()

uploads = parsed_batch.uploads
parsed_records = parsed_batch.parsed_records
parse_errors = parsed_batch.parse_errors

st.write(
    {
        "uploaded_files": len(uploads),
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

trim_scope = cast(
    TrimScope,
    st.radio(
        "Shared trim application mode",
        options=["all", "selected"],
        format_func=lambda value: (
            "All parsed records" if value == "all" else "Selected-record overrides"
        ),
        index=0 if trim_state.trim_scope == "all" else 1,
        horizontal=True,
        key=_BATCH_TRIM_SCOPE_WIDGET_KEY,
    ),
)
if trim_scope != trim_state.trim_scope:
    trim_state = BatchTrimState(
        trim_scope=trim_scope,
        global_trim_config=trim_state.global_trim_config,
        trim_configs_by_record=dict(trim_state.trim_configs_by_record),
    )
    from abi_sauce.viewer_state import set_batch_trim_state

    set_batch_trim_state(st.session_state, trim_state)

record_annotations = build_record_annotations(
    parsed_records.keys(),
    trim_state.trim_configs_by_record,
)

current_global_trim_config = trim_state.global_trim_config or TrimConfig()
should_refresh_trim_form = (
    st.session_state.get(_BATCH_TRIM_FORM_SCOPE_KEY) != trim_scope
    or st.session_state.get(_BATCH_TRIM_FORM_GLOBAL_CONFIG_KEY)
    != current_global_trim_config
)
if should_refresh_trim_form:
    _load_trim_config_into_form(current_global_trim_config)
    st.session_state[_BATCH_TRIM_FORM_SCOPE_KEY] = trim_scope
    st.session_state[_BATCH_TRIM_FORM_GLOBAL_CONFIG_KEY] = current_global_trim_config

st.caption(
    "The form below edits the shared all-record trim config. Per-sample overrides are edited in Sample Viewer."
)
if trim_scope == "selected":
    st.info(
        "Selected-record override mode is active. Batch results below reflect the per-sample overrides from Sample Viewer. "
        "The batch trim form still updates the stored all-record config for later use."
    )

with st.form("batch_trim_form"):
    st.subheader("Batch-wide trim defaults")
    left_trim = st.number_input(
        "Left trim",
        min_value=0,
        step=1,
        key=_BATCH_TRIM_LEFT_KEY,
    )
    right_trim = st.number_input(
        "Right trim",
        min_value=0,
        step=1,
        key=_BATCH_TRIM_RIGHT_KEY,
    )
    min_length = st.number_input(
        "Minimum length",
        min_value=0,
        step=1,
        key=_BATCH_TRIM_MIN_LENGTH_KEY,
    )
    quality_trim_enabled = st.checkbox(
        "Enable Mott quality trimming",
        key=_BATCH_TRIM_QUALITY_ENABLED_KEY,
    )
    st.number_input(
        "Mott max acceptable error probability",
        min_value=0.0,
        max_value=1.0,
        step=0.0001,
        format="%.4f",
        key=_BATCH_TRIM_ERROR_CUTOFF_KEY,
        disabled=not quality_trim_enabled,
        help=(
            "Lower values trim more aggressively. Examples: "
            "Q20 = 0.01, Q30 = 0.001, Q40 = 0.0001."
        ),
    )

    st.form_submit_button("Save batch trim defaults", on_click=_apply_global_trim_form)

trim_state = get_batch_trim_state(st.session_state)
resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)
prepared_batch = apply_trim_configs(
    parsed_batch,
    default_trim_config=resolved_trim_inputs.default_trim_config,
    trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
)
batch_summary = prepared_batch.batch_summary
batch_summary_rows = [
    {
        **record_summary.to_row(),
        "custom_trim": record_annotations.custom_trim_flags_by_record.get(
            record_summary.source_filename,
            False,
        ),
    }
    for record_summary in batch_summary.records
]

summary_col_1, summary_col_2, summary_col_3, summary_col_4 = st.columns(4)
with summary_col_1:
    st.metric("Trimmed records", batch_summary.trimmed_records)
with summary_col_2:
    st.metric("Passing minimum length", batch_summary.records_passing_min_length)
with summary_col_3:
    st.metric("Failing minimum length", batch_summary.records_failing_min_length)
with summary_col_4:
    st.metric("FASTQ exportable", batch_summary.fastq_exportable_records)

st.write(
    {
        "trim_scope": trim_state.trim_scope,
        "custom_trim_overrides": record_annotations.overridden_count,
        "batch_default_quality_trim_enabled": current_global_trim_config.quality_trim_enabled,
        "batch_default_error_probability_cutoff": current_global_trim_config.error_probability_cutoff,
    }
)

st.subheader("Batch summary")
st.dataframe(batch_summary_rows, hide_index=True, width="stretch")

st.subheader("Batch download")
export_format = st.selectbox(
    "Format",
    options=["fasta", "fastq"],
    key="batch_export_format",
)
export_format = cast(ExportFormat, export_format)

concatenate_batch = st.checkbox(
    "Concatenate entries into a single file",
    key="concatenate_batch",
)
batch_filename_stem = st.text_input(
    "Output filename stem",
    key="batch_filename_stem",
)

exclude_failed_min_length = st.checkbox(
    "Exclude sequences failing minimum length",
    key="exclude_failed_min_length_from_export",
)

try:
    download_artifact = prepare_batch_download(
        prepared_batch,
        export_format=export_format,
        concatenate_batch=bool(concatenate_batch),
        filename_stem=batch_filename_stem,
        require_min_length=bool(exclude_failed_min_length),
    )
except ExportError as exc:
    st.warning(str(exc))
else:
    st.write(
        {
            "exportable_records": len(download_artifact.eligible_records),
            "excluded_records": len(download_artifact.ineligible_reasons),
        }
    )

    if export_format == "fastq" and download_artifact.ineligible_reasons:
        st.info("FASTQ export will include only FASTQ-eligible records.")

    if download_artifact.ineligible_reasons:
        with st.expander("Excluded from current export"):
            for filename, reasons in download_artifact.ineligible_reasons:
                st.warning(f"{filename}: {', '.join(reasons)}")

    if not download_artifact.is_downloadable:
        st.warning("No trimmed records are eligible for this export selection.")
    else:
        st.download_button(
            label="Download trimmed batch",
            data=download_artifact.data,
            file_name=download_artifact.filename,
            mime=download_artifact.mime,
        )
