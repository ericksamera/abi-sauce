from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.chromatogram import build_chromatogram_view
from abi_sauce.chromatogram_figure import build_chromatogram_figure
from abi_sauce.services.batch import apply_trim_configs
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.trimming import TrimConfig

st.set_page_config(page_title="Chromatogram Preview", layout="wide")
st.title("Chromatogram Viewer")
st.caption("Inspect chromatograms from the shared active batch loaded in the sidebar.")

parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info("Upload one or more .ab1 files from the sidebar to inspect chromatograms.")
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

selected_record_name = cast(
    str,
    st.selectbox(
        "Selected record",
        options=list(parsed_records.keys()),
        format_func=lambda filename: (
            f"{filename} — "
            f"{parsed_records[filename].name or parsed_records[filename].record_id}"
        ),
    ),
)
record = parsed_records[selected_record_name]

st.subheader("Preview trim controls")
preview_scope = cast(
    str,
    st.radio(
        "Apply trim preview to",
        options=["selected", "all"],
        format_func=lambda value: (
            "Selected record only" if value == "selected" else "All parsed records"
        ),
        horizontal=True,
    ),
)

trim_left_col, trim_right_col, trim_min_length_col = st.columns(3)
with trim_left_col:
    left_trim = cast(
        int,
        st.number_input(
            "Left trim",
            min_value=0,
            step=1,
            value=0,
        ),
    )
with trim_right_col:
    right_trim = cast(
        int,
        st.number_input(
            "Right trim",
            min_value=0,
            step=1,
            value=0,
        ),
    )
with trim_min_length_col:
    min_length = cast(
        int,
        st.number_input(
            "Minimum length",
            min_value=0,
            step=1,
            value=1,
        ),
    )

quality_trim_enabled = cast(
    bool,
    st.checkbox(
        "Enable Mott quality trimming",
        value=False,
    ),
)
error_probability_cutoff = cast(
    float,
    st.number_input(
        "Mott max acceptable error probability",
        min_value=0.0,
        max_value=1.0,
        step=0.0001,
        format="%.4f",
        value=0.0100,
        disabled=not quality_trim_enabled,
        help=(
            "Lower values trim more aggressively. Examples: "
            "Q20 = 0.01, Q30 = 0.001, Q40 = 0.0001."
        ),
    ),
)

preview_trim_config = TrimConfig(
    left_trim=left_trim,
    right_trim=right_trim,
    min_length=min_length,
    quality_trim_enabled=quality_trim_enabled,
    error_probability_cutoff=error_probability_cutoff,
)

prepared_batch = apply_trim_configs(
    parsed_batch,
    default_trim_config=preview_trim_config if preview_scope == "all" else None,
    trim_configs_by_name=(
        {selected_record_name: preview_trim_config}
        if preview_scope == "selected"
        else None
    ),
)
trim_result = prepared_batch.trim_results[selected_record_name]
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
        "preview_scope": preview_scope,
        "quality_trim_enabled": quality_trim_enabled,
        "error_probability_cutoff": error_probability_cutoff,
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

with st.expander("Batch summary"):
    st.dataframe(
        prepared_batch.batch_summary.table_rows(),
        hide_index=True,
        width="stretch",
    )
