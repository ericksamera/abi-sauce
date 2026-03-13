from urllib.parse import quote

import streamlit as st

from abi_sauce.demo_batch import can_load_demo_sample, load_demo_sample
from abi_sauce.services.batch import apply_trim_configs, prepare_batch_download
from abi_sauce.trim_state import resolve_batch_trim_inputs
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import get_batch_trim_state


def _build_blast_url() -> str | None:
    parsed_batch = get_active_parsed_batch(st.session_state)
    if parsed_batch is None or not parsed_batch.parsed_records:
        return None

    trim_state = get_batch_trim_state(st.session_state)
    resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)
    prepared_batch = apply_trim_configs(
        parsed_batch,
        default_trim_config=resolved_trim_inputs.default_trim_config,
        trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
    )

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fasta",
        concatenate_batch=True,
        filename_stem="abi-sauce-batch",
        require_min_length=False,  # flip to True if you want to drop min-length failures
    )
    if not artifact.is_downloadable or not isinstance(artifact.data, str):
        return None

    query_text = artifact.data.strip()
    return (
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        "?PROGRAM=blastn"
        "&PAGE_TYPE=BlastSearch"
        "&LINK_LOC=blasthome"
        f"&QUERY={quote(query_text, safe='')}"
    )

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

st.subheader("BLAST")

blast_url = _build_blast_url()
if blast_url is None:
    st.caption("No trimmed sequences available for BLAST.")
elif len(blast_url) > 6000:
    st.caption(
        "Batch is too large for a reliable BLAST URL; export FASTA and paste/upload it in BLAST instead."
    )
else:
    st.link_button(
        "Open all trimmed sequences in blastn",
        blast_url,
        type="primary",
    )