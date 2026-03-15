from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.chromatogram import ChromatogramView, build_chromatogram_view
from abi_sauce.chromatogram_figure import build_chromatogram_figure
from abi_sauce.reference_alignment import (
    StrandPolicy,
    align_trimmed_read_to_reference,
    alignment_events_to_rows,
    format_alignment_block,
)
from abi_sauce.streamlit_cache import prepare_batch_for_trim_state
from abi_sauce.trim_state import build_record_annotations
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import (
    get_batch_trim_state,
    set_selected_record_name,
    sync_viewer_session_state,
)

_REFERENCE_INPUT_KEY = "reference_alignment.reference_input"
_REFERENCE_UPLOAD_KEY = "reference_alignment.reference_upload"
_STRAND_POLICY_KEY = "reference_alignment.strand_policy"
_SELECTED_RECORD_WIDGET_KEY = "reference_alignment.selected_record_widget"
_SELECTED_EVENT_WIDGET_KEY = "reference_alignment.selected_event_widget"

st.set_page_config(page_title="Reference Alignment", layout="wide")
st.title("Reference Alignment")
st.caption(
    "Align the currently trimmed read against a pasted reference sequence "
    "and inspect mismatches or indels against the chromatogram."
)


def _resolve_reference_text() -> str:
    pasted_reference = str(st.session_state.get(_REFERENCE_INPUT_KEY, "") or "").strip()
    if pasted_reference:
        return pasted_reference

    uploaded_reference = st.session_state.get(_REFERENCE_UPLOAD_KEY)
    if uploaded_reference is None:
        return ""

    try:
        return uploaded_reference.getvalue().decode("utf-8", errors="ignore")
    except AttributeError:
        return ""


def _average_base_spacing(view: ChromatogramView) -> float:
    steps = [
        float(right.position - left.position)
        for left, right in zip(view.base_calls, view.base_calls[1:])
        if right.position > left.position
    ]
    if steps:
        return sum(steps) / len(steps)
    return 10.0


def _centered_x_range(
    view: ChromatogramView,
    *,
    center: float,
    visible_bases: int = 25,
) -> list[float]:
    half_window = _average_base_spacing(view) * visible_bases
    left = max(0.0, center - half_window)
    right = min(float(max(view.trace_length - 1, 0)), center + half_window)
    if right <= left:
        right = min(float(max(view.trace_length - 1, 0)), left + max(half_window, 1.0))
    return [left, right]


def _event_option_label(event_row: dict[str, object]) -> str:
    event_type = str(event_row["type"])
    ref_pos = event_row["ref_pos"]
    query_pos = event_row["query_pos"]
    ref_base = str(event_row["ref_base"])
    query_base = str(event_row["query_base"])
    qscore = event_row["qscore"]
    qscore_label = "NA" if qscore is None else str(qscore)
    return (
        f"{event_type} | ref {ref_pos} | query {query_pos} | "
        f"{ref_base}>{query_base} | Q={qscore_label}"
    )


parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info("Load one or more .ab1 files from the workspace sidebar to begin.")
    st.stop()

parsed_records = parsed_batch.parsed_records
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
        "Sample",
        options=record_names,
        index=record_names.index(selected_record_name),
        format_func=lambda filename: record_annotations.display_labels_by_record[
            filename
        ],
        key=_SELECTED_RECORD_WIDGET_KEY,
    ),
)
set_selected_record_name(st.session_state, selected_record_name)

prepared_batch = prepare_batch_for_trim_state(
    parsed_batch,
    get_batch_trim_state(st.session_state),
)
raw_record = parsed_records[selected_record_name]
trim_result = prepared_batch.trim_results[selected_record_name]

if trim_result.trimmed_length <= 0:
    st.warning(
        "The selected sample currently trims down to an empty sequence. "
        "Relax the trim settings before aligning it to a reference."
    )
    st.stop()

reference_col, controls_col = st.columns([3, 1])

with reference_col:
    st.text_area(
        "Reference sequence",
        key=_REFERENCE_INPUT_KEY,
        height=180,
        placeholder="Paste a FASTA record or plain sequence here.",
    )
    st.file_uploader(
        "Or upload a FASTA/text reference",
        type=["fa", "fasta", "txt"],
        key=_REFERENCE_UPLOAD_KEY,
        accept_multiple_files=False,
    )

with controls_col:
    strand_policy = cast(
        StrandPolicy,
        st.selectbox(
            "Strand policy",
            options=["auto", "forward", "reverse-complement"],
            key=_STRAND_POLICY_KEY,
        ),
    )
    st.caption(
        "This prototype uses semi-global style scoring: internal "
        "mismatches and gaps are penalized, while terminal overhangs "
        "are tolerated."
    )
    st.caption(
        "Forward/reverse-complement strand choices are interpreted "
        "relative to the sample's displayed orientation: "
        f"{raw_record.orientation.replace('_', '-')}."
    )
    st.write(
        {
            "trimmed_query_length": trim_result.trimmed_length,
            "original_query_length": trim_result.original_length,
        }
    )

reference_text = _resolve_reference_text()
if not reference_text.strip():
    st.info("Paste or upload a reference sequence to run an alignment.")
    st.stop()

try:
    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=reference_text,
        strand_policy=strand_policy,
    )
except ValueError as exc:
    st.error(str(exc))
    st.stop()

metric_col_1, metric_col_2, metric_col_3, metric_col_4, metric_col_5, metric_col_6 = (
    st.columns(6)
)
with metric_col_1:
    st.metric("Strand", alignment_result.strand)
with metric_col_2:
    st.metric("Score", f"{alignment_result.score:.1f}")
with metric_col_3:
    st.metric("Identity", f"{alignment_result.percent_identity:.1f}%")
with metric_col_4:
    st.metric("Mismatches", alignment_result.mismatch_count)
with metric_col_5:
    st.metric("Insertions", alignment_result.insertion_count)
with metric_col_6:
    st.metric("Deletions", alignment_result.deletion_count)

st.caption(
    f"Reference: {alignment_result.reference_name} | "
    f"Reference span: {alignment_result.reference_start}–"
    f"{alignment_result.reference_end} | "
    f"Query span: {alignment_result.query_start}–{alignment_result.query_end}"
)

st.subheader("Gapped alignment")
st.code(format_alignment_block(alignment_result), wrap_lines=False)

event_rows = alignment_events_to_rows(alignment_result, include_matches=False)
selected_event_row: dict[str, object] | None = None

st.subheader("Alignment events")
if event_rows:
    st.caption(
        "Positions are 1-based within the aligned reference and the "
        "selected query strand relative to the sample's displayed orientation."
    )
    st.dataframe(event_rows, hide_index=True, width="stretch")
    selected_event_index = cast(
        int,
        st.selectbox(
            "Center chromatogram on event",
            options=list(range(len(event_rows))),
            format_func=lambda index: _event_option_label(event_rows[index]),
            key=_SELECTED_EVENT_WIDGET_KEY,
        ),
    )
    selected_event_row = event_rows[selected_event_index]
else:
    st.success(
        "No mismatch or indel events were detected in the best-scoring alignment."
    )

chromatogram_view = build_chromatogram_view(raw_record, trim_result)
if not chromatogram_view.is_renderable:
    st.warning("Selected record does not have enough trace data to render a chart.")
    st.stop()

theme_type = str(getattr(getattr(st.context, "theme", None), "type", "light"))
figure = build_chromatogram_figure(
    chromatogram_view,
    theme_type=theme_type,
)
trace_x = None if selected_event_row is None else selected_event_row.get("trace_x")
if isinstance(trace_x, (int, float)):
    figure.update_xaxes(
        range=_centered_x_range(
            chromatogram_view,
            center=float(trace_x),
        )
    )
figure.update_layout(
    height=500,
    margin={"l": 24, "r": 24, "t": 24, "b": 24},
)

st.subheader("Chromatogram")
if selected_event_row is not None:
    if isinstance(trace_x, (int, float)):
        st.caption(f"Centered on trace sample {trace_x} for the selected event.")
    else:
        st.caption(
            "The selected event is a deletion without a direct query peak, "
            "so the chromatogram remains at the default view."
        )
st.plotly_chart(
    figure,
    width="stretch",
    config={"scrollZoom": False},
)
