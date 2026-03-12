from __future__ import annotations

import streamlit as st

from abi_sauce.exceptions import AbiParseError
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

with st.form("trim_form"):
    st.subheader("Trim")
    left_trim = st.number_input("Left trim", min_value=0, value=0, step=1)
    right_trim = st.number_input("Right trim", min_value=0, value=0, step=1)
    min_length = st.number_input("Minimum length", min_value=0, value=1, step=1)
    apply_trim = st.form_submit_button("Apply trim")

if not apply_trim:
    st.subheader("Raw sequence preview")
    st.code(record.sequence[:500] or "<empty>")
    st.stop()

trim_result = trim_sequence_record(
    record,
    TrimConfig(
        left_trim=int(left_trim),
        right_trim=int(right_trim),
        min_length=int(min_length),
    ),
)

col1, col2 = st.columns(2)
with col1:
    st.metric("Original length", trim_result.original_length)
with col2:
    st.metric("Trimmed length", trim_result.trimmed_length)

st.write(
    {
        "bases_removed_left": trim_result.bases_removed_left,
        "bases_removed_right": trim_result.bases_removed_right,
        "passed_min_length": trim_result.passed_min_length,
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
