from __future__ import annotations

import streamlit as st

from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.demo_batch import can_load_demo_sample, load_demo_sample

st.set_page_config(page_title="ABI Sauce", layout="wide")
st.title("ABI Sauce")
st.caption(
    "Load a shared batch once, inspect one sample in the Sample Viewer, and make batch-wide changes in the Batch Viewer."
)

parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info(
        "Use the sidebar to load one or more ABI trace files into the current workspace."
    )
    st.subheader("Suggested workflow")
    st.markdown(
        "\n".join(
            [
                "1. Load ABI trace files from the sidebar.",
                "2. Open **Sample Viewer** to inspect one selected sample and its chromatogram.",
                "3. Open **Batch Viewer** to apply batch-wide trim settings and export processed records.",
            ]
        )
    )

    demo_available = can_load_demo_sample()

    if st.button(
        "Load demo sample",
        disabled=not demo_available,
    ):
        load_demo_sample(st.session_state)
        st.rerun()

    if not demo_available:
        st.caption("Demo sample unavailable: tests/fixtures/example.ab1 was not found.")
    st.stop()

uploaded_count = len(parsed_batch.uploads)
parsed_count = len(parsed_batch.parsed_records)
failed_count = len(parsed_batch.parse_errors)

metric_col_1, metric_col_2, metric_col_3 = st.columns(3)
with metric_col_1:
    st.metric("Uploaded files", uploaded_count)
with metric_col_2:
    st.metric("Parsed records", parsed_count)
with metric_col_3:
    st.metric("Parse failures", failed_count)

viewer_col, batch_col = st.columns(2)
with viewer_col:
    st.subheader("Sample Viewer")
    st.write(
        "Use the Sample Viewer to inspect one selected sample, visualize chromatograms, and edit per-sample trim behavior."
    )
with batch_col:
    st.subheader("Batch Viewer")
    st.write(
        "Use the Batch Viewer to manage shared batch-wide trim settings, review all parsed records, and prepare exports."
    )

if parsed_batch.parsed_records:
    st.subheader("Parsed records")
    st.dataframe(
        [
            {
                "filename": filename,
                "record_id": record.record_id,
                "name": record.name,
                "sequence_length": len(record.sequence),
                "has_trace_data": record.trace_data is not None,
                "quality_count": (
                    0 if record.qualities is None else len(record.qualities)
                ),
            }
            for filename, record in parsed_batch.parsed_records.items()
        ],
        hide_index=True,
        width="stretch",
    )

if parsed_batch.parse_errors:
    with st.expander("Files that failed to parse"):
        for filename, message in parsed_batch.parse_errors.items():
            st.error(f"{filename}: {message}")
