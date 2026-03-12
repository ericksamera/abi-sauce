from __future__ import annotations

import streamlit as st

from abi_sauce.exceptions import AbiParseError
from abi_sauce.export import ExportError, to_fasta, to_fastq
from abi_sauce.models import SequenceUpload
from abi_sauce.parsers.abi import parse_ab1_upload
from abi_sauce.trimming import TrimConfig, trim_sequence_record

st.set_page_config(page_title="ABI Sauce", layout="wide")
st.title("ABI Sauce")

uploaded_file = st.file_uploader(
    "Upload one ABI trace file",
    type=["ab1", "abi"],
)

if uploaded_file is None:
    st.info("Upload one .ab1 file to test the parser.")
    st.stop()

upload = SequenceUpload(
    filename=uploaded_file.name,
    content=uploaded_file.getvalue(),
)

st.write(
    {
        "filename": upload.filename,
        "suffix": upload.suffix,
        "size_bytes": upload.size_bytes,
    }
)

try:
    record = parse_ab1_upload(upload)
except AbiParseError as exc:
    st.error(str(exc))
    st.stop()
except Exception as exc:
    st.error("Parsing failed.")
    st.exception(exc)
    st.stop()

st.success("Parsed successfully.")
st.write(
    {
        "record_id": record.record_id,
        "name": record.name,
        "source_format": record.source_format,
        "sequence_length": len(record.sequence),
        "quality_count": 0 if record.qualities is None else len(record.qualities),
        "has_trace_data": record.trace_data is not None,
    }
)

# Reset applied trim when a new upload appears.
upload_signature = (upload.filename, upload.size_bytes)
if st.session_state.get("upload_signature") != upload_signature:
    st.session_state.upload_signature = upload_signature
    st.session_state.applied_trim_config = None
    st.session_state.trim_quality_enabled = False
    st.session_state.trim_error_probability_cutoff = 0.01

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

trim_result = trim_sequence_record(record, applied_trim_config)

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

st.subheader("Download")
export_format = st.selectbox(
    "Format",
    options=["fasta", "fastq"],
    key="export_format",
)

try:
    if export_format == "fasta":
        export_text = to_fasta(trim_result.record)
        export_filename = (
            f"{trim_result.record.name or trim_result.record.record_id}.fasta"
        )
    else:
        export_text = to_fastq(trim_result.record)
        export_filename = (
            f"{trim_result.record.name or trim_result.record.record_id}.fastq"
        )
except ExportError as exc:
    st.warning(str(exc))
else:
    st.download_button(
        label="Download trimmed sequence",
        data=export_text,
        file_name=export_filename,
        mime="text/plain",
    )
