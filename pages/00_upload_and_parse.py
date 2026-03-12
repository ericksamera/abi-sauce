from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.batch import ExportFormat
from abi_sauce.services.batch import (
    apply_trim_config,
    parse_uploaded_batch,
    prepare_batch_download,
)
from abi_sauce.exceptions import ExportError
from abi_sauce.trimming import TrimConfig

st.set_page_config(page_title="ABI Sauce", layout="wide")
st.title("ABI Sauce")

uploaded_files = st.file_uploader(
    "Upload ABI trace files",
    type=["ab1", "abi"],
    accept_multiple_files=True,
)

if not uploaded_files:
    st.info("Upload one or more .ab1 files to test the parser.")
    st.stop()

parsed_batch = parse_uploaded_batch(uploaded_files)
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

if st.session_state.get("upload_signature") != parsed_batch.signature:
    st.session_state.upload_signature = parsed_batch.signature
    st.session_state.applied_trim_config = None
    st.session_state.selected_record_name = next(iter(parsed_records))
    st.session_state.trim_quality_enabled = False
    st.session_state.trim_error_probability_cutoff = 0.01
    st.session_state.batch_export_format = "fasta"
    st.session_state.concatenate_batch = True
    st.session_state.batch_filename_stem = "abi-sauce-batch"
    st.session_state.exclude_failed_min_length_from_export = True

selected_record_name = st.selectbox(
    "Selected record",
    options=list(parsed_records.keys()),
    key="selected_record_name",
)
record = parsed_records[selected_record_name]

st.success("Parsed successfully.")
st.write(
    {
        "selected_record": selected_record_name,
        "record_id": record.record_id,
        "name": record.name,
        "source_format": record.source_format,
        "sequence_length": len(record.sequence),
        "quality_count": 0 if record.qualities is None else len(record.qualities),
        "has_trace_data": record.trace_data is not None,
    }
)

with st.form("trim_form"):
    st.subheader("Trim")
    left_trim = st.number_input(
        "Left trim",
        min_value=0,
        value=0,
        step=1,
        key="trim_left",
    )
    right_trim = st.number_input(
        "Right trim",
        min_value=0,
        value=0,
        step=1,
        key="trim_right",
    )
    min_length = st.number_input(
        "Minimum length",
        min_value=0,
        value=1,
        step=1,
        key="trim_min_length",
    )
    quality_trim_enabled = st.checkbox(
        "Enable Mott quality trimming",
        key="trim_quality_enabled",
    )
    error_probability_cutoff = st.number_input(
        "Mott max acceptable error probability",
        min_value=0.0,
        max_value=1.0,
        step=0.0001,
        format="%.4f",
        key="trim_error_probability_cutoff",
        disabled=not quality_trim_enabled,
        help=(
            "Lower values trim more aggressively. Examples: "
            "Q20 = 0.01, Q30 = 0.001, Q40 = 0.0001."
        ),
    )
    submitted = st.form_submit_button("Apply trim")

if submitted:
    st.session_state.applied_trim_config = TrimConfig(
        left_trim=int(left_trim),
        right_trim=int(right_trim),
        min_length=int(min_length),
        quality_trim_enabled=bool(quality_trim_enabled),
        error_probability_cutoff=float(error_probability_cutoff),
    )

applied_trim_config = st.session_state.get("applied_trim_config")

if applied_trim_config is None:
    st.subheader("Raw sequence preview")
    st.code(record.sequence[:500] or "<empty>")
    st.stop()

prepared_batch = apply_trim_config(parsed_batch, applied_trim_config)
trim_results = prepared_batch.trim_results
batch_summary = prepared_batch.batch_summary
trim_result = trim_results[selected_record_name]

st.subheader("Batch summary")
st.write(
    {
        "trimmed_records": batch_summary.trimmed_records,
        "records_passing_min_length": batch_summary.records_passing_min_length,
        "records_failing_min_length": batch_summary.records_failing_min_length,
        "fastq_exportable_records": batch_summary.fastq_exportable_records,
    }
)
st.dataframe(batch_summary.table_rows(), hide_index=True, width="stretch")

st.subheader("Selected record detail")

col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Original length", trim_result.original_length)
with col2:
    st.metric("Trimmed length", trim_result.trimmed_length)
with col3:
    st.metric(
        "Quality trim left",
        trim_result.quality_bases_removed_left,
    )
with col4:
    st.metric(
        "Quality trim right",
        trim_result.quality_bases_removed_right,
    )

st.write(
    {
        "bases_removed_left": trim_result.bases_removed_left,
        "bases_removed_right": trim_result.bases_removed_right,
        "passed_min_length": trim_result.passed_min_length,
        "quality_trim_enabled": applied_trim_config.quality_trim_enabled,
        "error_probability_cutoff": applied_trim_config.error_probability_cutoff,
        "batch_records": len(trim_results),
    }
)

if not trim_result.passed_min_length:
    st.warning("Trimmed sequence did not meet the minimum length.")

raw_col, trimmed_col = st.columns(2)
with raw_col:
    st.subheader("Raw sequence")
    st.code(record.sequence[:500] or "<empty>")

with trimmed_col:
    st.subheader("Trimmed sequence")
    st.code(trim_result.record.sequence[:500] or "<empty>")

st.subheader("Batch download")
export_format = st.selectbox(
    "Format",
    options=["fasta", "fastq"],
    key="batch_export_format",
)
export_format = cast(ExportFormat, export_format)

concatenate_batch = st.checkbox(
    "Concatenate entries into a single file",
    value=True,
    key="concatenate_batch",
)
batch_filename_stem = st.text_input(
    "Output filename stem",
    value="abi-sauce-batch",
    key="batch_filename_stem",
)

exclude_failed_min_length = st.checkbox(
    "Exclude sequences failing minimum length",
    value=True,
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
