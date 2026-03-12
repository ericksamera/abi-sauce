from __future__ import annotations

import streamlit as st

from abi_sauce.exceptions import AbiParseError
from abi_sauce.models import SequenceUpload
from abi_sauce.parsers.abi import parse_ab1_upload

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
