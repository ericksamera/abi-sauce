import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from abi_sauce.demo_batch import can_load_demo_sample, load_demo_sample
from abi_sauce.models import SequenceRecord
from abi_sauce.services.batch import apply_trim_configs
from abi_sauce.trim_state import resolve_batch_trim_inputs
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import get_batch_trim_state

_QC_METRICS: tuple[tuple[str, str], ...] = (
    ("trace_score", "Trace score"),
    ("pup_score", "PUP score"),
    ("crl_score", "CRL"),
)
_LENGTH_SERIES: tuple[tuple[str, str], ...] = (
    ("original_length", "Before trim"),
    ("trimmed_length", "After trim"),
)


def _annotation_number(record: SequenceRecord, key: str) -> float | None:
    value = record.annotations.get(key)
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _build_record_rows(
    parsed_records: dict[str, SequenceRecord],
) -> list[dict[str, object]]:
    return [
        {
            "filename": filename,
            "record_id": record.record_id,
            "name": record.name,
            "sequence_length": len(record.sequence),
            "orientation": record.orientation,
            "has_trace_data": record.trace_data is not None,
            "quality_count": 0 if record.qualities is None else len(record.qualities),
            "trace_score": _annotation_number(record, "trace_score"),
            "pup_score": _annotation_number(record, "pup_score"),
            "crl_score": _annotation_number(record, "crl_score"),
        }
        for filename, record in parsed_records.items()
    ]


def _build_length_rows(
    parsed_records: dict[str, SequenceRecord],
    trimmed_lengths_by_filename: dict[str, int],
) -> list[dict[str, object]]:
    return [
        {
            "filename": filename,
            "record_id": record.record_id,
            "name": record.name,
            "original_length": len(record.sequence),
            "trimmed_length": trimmed_lengths_by_filename.get(
                filename,
                len(record.sequence),
            ),
        }
        for filename, record in parsed_records.items()
    ]


def _sample_label(row: dict[str, object]) -> str:
    for key in ("name", "filename", "record_id"):
        value = row.get(key)
        if isinstance(value, str) and value.strip():
            return value
    return "unknown sample"


def _metric_values(record_rows: list[dict[str, object]], key: str) -> list[float]:
    values: list[float] = []
    for row in record_rows:
        value = row.get(key)
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)):
            values.append(float(value))
    return values


def _metric_points(
    record_rows: list[dict[str, object]],
    key: str,
) -> tuple[list[float], list[str]]:
    values: list[float] = []
    sample_labels: list[str] = []
    for row in record_rows:
        value = row.get(key)
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)):
            values.append(float(value))
            sample_labels.append(_sample_label(row))
    return values, sample_labels


def _build_qc_distribution_figure(record_rows: list[dict[str, object]]) -> go.Figure:
    figure = make_subplots(
        rows=2,
        cols=len(_QC_METRICS),
        subplot_titles=[label for _, label in _QC_METRICS],
        shared_xaxes="columns",
        row_heights=[0.22, 0.78],
        horizontal_spacing=0.08,
        vertical_spacing=0.05,
    )

    for column_index, (key, label) in enumerate(_QC_METRICS, start=1):
        values, sample_labels = _metric_points(record_rows, key)
        if not values:
            continue

        figure.add_trace(
            go.Box(
                x=values,
                name=label,
                showlegend=False,
                orientation="h",
                boxpoints="outliers",
                marker={"size": 5},
                customdata=sample_labels,
                hovertemplate=(
                    "sample=%{customdata}<br>" f"{label}=%{{x}}" "<extra></extra>"
                ),
            ),
            row=1,
            col=column_index,
        )
        figure.add_trace(
            go.Histogram(
                x=values,
                name=label,
                showlegend=False,
                nbinsx=max(5, min(20, len(values))),
                hovertemplate=f"{label}=%{{x}}<br>count=%{{y}}<extra></extra>",
            ),
            row=2,
            col=column_index,
        )
        figure.update_xaxes(title_text=label, row=2, col=column_index)
        figure.update_yaxes(visible=False, row=1, col=column_index)
        figure.update_yaxes(title_text="Count", row=2, col=column_index)

    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_layout(
        height=280,
        margin={"l": 16, "r": 16, "t": 44, "b": 16},
        bargap=0.1,
    )
    return figure


def _build_length_distribution_figure(
    length_rows: list[dict[str, object]],
) -> go.Figure:
    figure = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.22, 0.78],
        vertical_spacing=0.05,
    )
    nbins = max(10, min(40, len(length_rows)))

    for key, label in _LENGTH_SERIES:
        values, sample_labels = _metric_points(length_rows, key)
        if not values:
            continue

        figure.add_trace(
            go.Scatter(
                x=values,
                y=[label] * len(values),
                mode="markers",
                name=label,
                showlegend=False,
                legendgroup=key,
                marker={
                    "symbol": "line-ns-open",
                    "size": 18,
                    "line": {"width": 1.5},
                },
                customdata=sample_labels,
                hovertemplate=(
                    "sample=%{customdata}<br>"
                    f"state={label}<br>"
                    "length=%{x}<extra></extra>"
                ),
            ),
            row=1,
            col=1,
        )
        figure.add_trace(
            go.Histogram(
                x=values,
                name=label,
                showlegend=True,
                legendgroup=key,
                bingroup="sequence_length_distribution",
                nbinsx=nbins,
                opacity=0.7,
                hovertemplate=(
                    f"state={label}<br>" "length=%{x}<br>" "count=%{y}<extra></extra>"
                ),
            ),
            row=2,
            col=1,
        )

    figure.update_xaxes(showticklabels=False, row=1, col=1)
    figure.update_yaxes(
        title_text="",
        showgrid=False,
        zeroline=False,
        row=1,
        col=1,
    )
    figure.update_xaxes(title_text="Sequence length", row=2, col=1)
    figure.update_yaxes(title_text="Count", row=2, col=1)
    figure.update_layout(
        height=300,
        margin={"l": 16, "r": 16, "t": 24, "b": 16},
        barmode="overlay",
        hovermode="closest",
    )
    return figure


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

trim_state = get_batch_trim_state(st.session_state)
resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)
prepared_batch = apply_trim_configs(
    parsed_batch,
    default_trim_config=resolved_trim_inputs.default_trim_config,
    trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
)


uploaded_count = len(parsed_batch.uploads)
parsed_count = len(parsed_batch.parsed_records)
failed_count = len(parsed_batch.parse_errors)
record_rows = _build_record_rows(parsed_batch.parsed_records)
length_rows = _build_length_rows(
    parsed_batch.parsed_records,
    {
        filename: trim_result.trimmed_length
        for filename, trim_result in prepared_batch.trim_results.items()
    },
)
qc_rows_present = any(_metric_values(record_rows, key) for key, _ in _QC_METRICS)

metric_col_1, metric_col_2, metric_col_3 = st.columns(3)
with metric_col_1:
    st.metric("Uploaded files", uploaded_count)
with metric_col_2:
    st.metric("Parsed records", parsed_count)
with metric_col_3:
    st.metric("Parse failures", failed_count)

if length_rows:
    st.subheader("Sequence length distribution")
    st.caption(
        "Original and trimmed lengths are overlaid so you can compare how the current batch trim settings shift the length distribution across all loaded samples."
    )
    st.plotly_chart(
        _build_length_distribution_figure(length_rows),
        width="stretch",
        config={"scrollZoom": False, "displayModeBar": False},
    )

if parsed_batch.parsed_records and qc_rows_present:
    st.subheader("ABI QC score distributions")
    st.caption(
        "Each panel uses a top marginal boxplot over a histogram so you can see both spread and shape for the currently loaded samples. Missing ABI QC values are excluded metric-by-metric."
    )
    st.plotly_chart(
        _build_qc_distribution_figure(record_rows),
        width="stretch",
        config={"scrollZoom": False, "displayModeBar": False},
    )

if parsed_batch.parsed_records:
    st.subheader("Parsed records")
    st.dataframe(
        record_rows,
        hide_index=True,
        width="stretch",
        column_config={
            "orientation": st.column_config.TextColumn("Orientation"),
            "trace_score": st.column_config.NumberColumn(
                "Trace score",
                format="%.0f",
            ),
            "pup_score": st.column_config.NumberColumn(
                "PUP score",
                format="%.0f",
            ),
            "crl_score": st.column_config.NumberColumn(
                "CRL",
                format="%.0f",
            ),
        },
    )

if parsed_batch.parse_errors:
    with st.expander("Files that failed to parse"):
        for filename, message in parsed_batch.parse_errors.items():
            st.error(f"{filename}: {message}")
