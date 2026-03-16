from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Protocol, cast

import streamlit as st

from abi_sauce.alignment_state import (
    AlignmentDefinition,
    AlignmentEngineKind,
    create_alignment_definition,
    delete_alignment_definition,
    set_selected_alignment_id,
    suggest_alignment_name,
    sync_alignment_session_state,
    update_alignment_definition,
)
from abi_sauce.assembly_exports import format_assembly_alignment_fasta
from abi_sauce.assembly_presenters import (
    assembly_conflicts_to_rows,
    format_assembly_block,
)
from abi_sauce.assembly_trace_figure import build_assembly_trace_figure
from abi_sauce.assembly_types import AssemblyConfig, AssemblyResult
from abi_sauce.export import to_fasta
from abi_sauce.reference_alignment import normalize_reference
from abi_sauce.reference_alignment_exports import format_reference_alignment_fasta
from abi_sauce.reference_alignment_multi_exports import (
    consensus_record_from_reference_multi_result,
    format_reference_multi_alignment_fasta,
)
from abi_sauce.reference_alignment_multi_presenters import (
    format_reference_multi_alignment_block,
)
from abi_sauce.reference_alignment_multi_trace_figure import (
    build_reference_multi_alignment_trace_figure,
)
from abi_sauce.reference_alignment_presenters import format_alignment_block
from abi_sauce.reference_alignment_trace_figure import (
    build_reference_alignment_trace_figure,
)
from abi_sauce.reference_alignment_types import StrandPolicy
from abi_sauce.reference_library_state import (
    StoredReference,
    get_stored_reference,
    list_stored_references,
    store_reference,
)
from abi_sauce.services.alignment_compute import (
    ComputedAlignment,
    compute_saved_alignments,
)
from abi_sauce.streamlit_cache import (
    build_selected_assembly_trace_view,
    prepare_batch_for_trim_state,
)
from abi_sauce.trim_state import build_record_annotations
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import get_batch_trim_state, sync_viewer_session_state

_ALIGNMENT_SELECT_WIDGET_KEY = "alignments.selected_alignment_widget"
_SELECTED_EVENT_WIDGET_KEY = "alignments.selected_event_widget"

_NEW_NAME_WIDGET_KEY = "alignments.dialog.new.name"
_NEW_READS_WIDGET_KEY = "alignments.dialog.new.reads"
_NEW_ENGINE_WIDGET_KEY = "alignments.dialog.new.engine_kind"
_NEW_REFERENCE_SOURCE_WIDGET_KEY = "alignments.dialog.new.reference_source"
_NEW_EXISTING_REFERENCE_ID_WIDGET_KEY = "alignments.dialog.new.reference_existing_id"
_NEW_REFERENCE_INPUT_MODE_WIDGET_KEY = "alignments.dialog.new.reference_input_mode"
_NEW_REFERENCE_NAME_WIDGET_KEY = "alignments.dialog.new.reference_name"
_NEW_REFERENCE_TEXT_WIDGET_KEY = "alignments.dialog.new.reference_text"
_NEW_REFERENCE_UPLOAD_WIDGET_KEY = "alignments.dialog.new.reference_upload"
_NEW_STRAND_POLICY_WIDGET_KEY = "alignments.dialog.new.strand_policy"
_NEW_MIN_OVERLAP_WIDGET_KEY = "alignments.dialog.new.min_overlap"
_NEW_MIN_IDENTITY_WIDGET_KEY = "alignments.dialog.new.min_identity"
_NEW_QUALITY_MARGIN_WIDGET_KEY = "alignments.dialog.new.quality_margin"

_EDIT_NAME_WIDGET_KEY_PREFIX = "alignments.dialog.edit.name"
_EDIT_READS_WIDGET_KEY_PREFIX = "alignments.dialog.edit.reads"
_EDIT_ENGINE_WIDGET_KEY_PREFIX = "alignments.dialog.edit.engine_kind"
_EDIT_REFERENCE_SOURCE_WIDGET_KEY_PREFIX = "alignments.dialog.edit.reference_source"
_EDIT_EXISTING_REFERENCE_ID_WIDGET_KEY_PREFIX = (
    "alignments.dialog.edit.reference_existing_id"
)
_EDIT_REFERENCE_INPUT_MODE_WIDGET_KEY_PREFIX = (
    "alignments.dialog.edit.reference_input_mode"
)
_EDIT_REFERENCE_NAME_WIDGET_KEY_PREFIX = "alignments.dialog.edit.reference_name"
_EDIT_REFERENCE_TEXT_WIDGET_KEY_PREFIX = "alignments.dialog.edit.reference_text"
_EDIT_REFERENCE_UPLOAD_WIDGET_KEY_PREFIX = "alignments.dialog.edit.reference_upload"
_EDIT_STRAND_POLICY_WIDGET_KEY_PREFIX = "alignments.dialog.edit.strand_policy"
_EDIT_MIN_OVERLAP_WIDGET_KEY_PREFIX = "alignments.dialog.edit.min_overlap"
_EDIT_MIN_IDENTITY_WIDGET_KEY_PREFIX = "alignments.dialog.edit.min_identity"
_EDIT_QUALITY_MARGIN_WIDGET_KEY_PREFIX = "alignments.dialog.edit.quality_margin"

ReferenceSourceMode = Literal["existing", "new"]
NewReferenceInputMode = Literal["upload", "paste"]

_ALIGNED_TRACE_ROW_HEIGHT_PX = 84
_ALIGNED_TRACE_BASE_HEIGHT_PX = 72
_ALIGNED_TRACE_MIN_HEIGHT_PX = 240
_ALIGNED_TRACE_VISIBLE_COLUMNS = 150

st.set_page_config(page_title="Alignments", layout="wide")
st.title("Alignments")
st.caption(
    "Manage saved pairwise and reference-guided alignment definitions, "
    "recompute them against the current trimmed batch, and inspect the results."
)


def _engine_label(engine_kind: AlignmentEngineKind) -> str:
    if engine_kind == "pairwise":
        return "Pairwise assembly"
    if engine_kind == "reference_single":
        return "One read vs reference"
    return "Many reads vs reference"


def _alignment_member_count_valid(
    engine_kind: AlignmentEngineKind,
    *,
    selected_read_count: int,
) -> bool:
    if engine_kind == "pairwise":
        return selected_read_count == 2
    if engine_kind == "reference_single":
        return selected_read_count == 1
    return selected_read_count >= 2


def _reference_required(engine_kind: AlignmentEngineKind) -> bool:
    return engine_kind in {"reference_single", "reference_multi"}


def _centered_alignment_x_range(
    *,
    alignment_length: int,
    cell_width: float,
    center_column_index: int | None,
    visible_columns: int = _ALIGNED_TRACE_VISIBLE_COLUMNS,
) -> list[float]:
    if alignment_length <= 0:
        return [0.0, 1.0]

    full_left = 0.0
    full_right = float(alignment_length) * cell_width
    visible_width = min(float(visible_columns) * cell_width, full_right - full_left)
    if visible_width <= 0:
        return [full_left, full_right]

    if center_column_index is None:
        return [full_left, full_left + visible_width]

    center = (float(center_column_index) - 0.5) * cell_width
    half_width = visible_width / 2.0
    left = max(full_left, center - half_width)
    right = min(full_right, center + half_width)
    if right - left < visible_width:
        if left <= full_left:
            right = min(full_right, left + visible_width)
        elif right >= full_right:
            left = max(full_left, right - visible_width)
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


def _alignment_select_label(
    alignment_definitions_by_id: dict[str, AlignmentDefinition],
    alignment_id: str,
) -> str:
    definition = alignment_definitions_by_id[alignment_id]
    return definition.name


def _safe_filename_stem(value: str) -> str:
    collapsed = "_".join(value.split())
    cleaned = "".join(
        character
        for character in collapsed
        if character.isalnum() or character in {"-", "_", "."}
    ).strip("._-")
    return cleaned or "alignment"


def _default_dialog_name(
    *,
    source_filenames: tuple[str, ...],
    engine_kind: AlignmentEngineKind,
    reference_name: str | None,
    existing_names: tuple[str, ...],
    current_name: str | None = None,
) -> str:
    if current_name is not None and current_name.strip():
        return current_name
    return suggest_alignment_name(
        source_filenames,
        engine_kind=engine_kind,
        reference_name=reference_name,
        existing_names=existing_names,
    )


def _normalized_optional_text(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    stripped_value = value.strip()
    return stripped_value or None


def _selection_index(
    options: tuple[str, ...],
    *,
    current_value: str | None,
    default_value: str,
) -> int:
    resolved_value = (
        current_value
        if current_value in options
        else (default_value if default_value in options else options[0])
    )
    return options.index(resolved_value)


class _UploadedReferenceLike(Protocol):
    name: str

    def getvalue(self) -> bytes: ...


def _coerce_uploaded_reference(
    uploaded_reference: object,
) -> _UploadedReferenceLike | None:
    if uploaded_reference is None:
        return None
    if not hasattr(uploaded_reference, "getvalue"):
        return None
    if not hasattr(uploaded_reference, "name"):
        return None
    return cast(_UploadedReferenceLike, uploaded_reference)


def _uploaded_reference_text(uploaded_reference: object) -> str:
    uploaded_file = _coerce_uploaded_reference(uploaded_reference)
    if uploaded_file is None:
        return ""
    uploaded_bytes = uploaded_file.getvalue()
    if not isinstance(uploaded_bytes, bytes):
        return ""
    return uploaded_bytes.decode("utf-8", errors="ignore")


def _uploaded_reference_name(uploaded_reference: object) -> str | None:
    uploaded_file = _coerce_uploaded_reference(uploaded_reference)
    if uploaded_file is None:
        return None
    uploaded_name = uploaded_file.name
    if not isinstance(uploaded_name, str) or not uploaded_name.strip():
        return None
    return Path(uploaded_name).stem or None


def _raw_reference_has_header(reference_text: str) -> bool:
    for line in reference_text.splitlines():
        if line.strip():
            return line.lstrip().startswith(">")
    return False


def _canonical_reference_input(
    *,
    reference_name: str | None,
    reference_text: str,
    fallback_name: str | None = None,
) -> tuple[str, str]:
    parsed_name, sequence = normalize_reference(reference_text)
    explicit_name = _normalized_optional_text(reference_name)
    resolved_name = (
        explicit_name
        or (
            parsed_name
            if _raw_reference_has_header(reference_text)
            else _normalized_optional_text(fallback_name)
        )
        or parsed_name
    )
    return resolved_name, f">{resolved_name}\n{sequence}\n"


def _reference_display_label(stored_reference: StoredReference) -> str:
    _reference_name, reference_sequence = normalize_reference(
        stored_reference.reference_text
    )
    return f"{stored_reference.name} ({len(reference_sequence)} bp)"


def _matching_stored_reference_id(
    stored_references: tuple[StoredReference, ...],
    *,
    reference_name: str | None,
    reference_text: str | None,
) -> str | None:
    normalized_reference_text = _normalized_optional_text(reference_text)
    if normalized_reference_text is None:
        return None
    try:
        _resolved_name, canonical_reference_text = _canonical_reference_input(
            reference_name=reference_name,
            reference_text=normalized_reference_text,
        )
    except ValueError:
        return None

    for stored_reference in stored_references:
        if stored_reference.reference_text == canonical_reference_text:
            return stored_reference.reference_id
    return None


@dataclass(frozen=True, slots=True)
class _ResolvedReferenceSelection:
    source_mode: ReferenceSourceMode
    reference_name: str | None
    reference_text: str | None
    is_valid: bool


def _render_reference_selection(
    *,
    engine_kind: AlignmentEngineKind,
    stored_references: tuple[StoredReference, ...],
    reference_source_key: str,
    existing_reference_id_key: str,
    reference_input_mode_key: str,
    reference_name_key: str,
    reference_text_key: str,
    reference_upload_key: str,
    default_reference_name: str | None = None,
    default_reference_text: str | None = None,
) -> _ResolvedReferenceSelection:
    if not _reference_required(engine_kind):
        return _ResolvedReferenceSelection(
            source_mode="new",
            reference_name=None,
            reference_text=None,
            is_valid=True,
        )

    matched_reference_id = _matching_stored_reference_id(
        stored_references,
        reference_name=default_reference_name,
        reference_text=default_reference_text,
    )
    source_options: tuple[ReferenceSourceMode, ...]
    if stored_references:
        source_options = ("existing", "new")
    else:
        source_options = ("new",)

    default_source_mode: ReferenceSourceMode = (
        "existing"
        if stored_references
        and (
            matched_reference_id is not None
            or _normalized_optional_text(default_reference_text) is None
        )
        else "new"
    )
    reference_source = cast(
        ReferenceSourceMode,
        st.selectbox(
            "Reference source",
            options=list(source_options),
            index=_selection_index(
                source_options,
                current_value=cast(
                    str | None,
                    st.session_state.get(reference_source_key),
                ),
                default_value=default_source_mode,
            ),
            format_func=lambda value: (
                "Use existing reference"
                if value == "existing"
                else "Upload/paste new reference"
            ),
            key=reference_source_key,
        ),
    )

    if reference_source == "existing":
        reference_ids = tuple(
            stored_reference.reference_id for stored_reference in stored_references
        )
        selected_reference_id = cast(
            str,
            st.selectbox(
                "Existing reference",
                options=list(reference_ids),
                index=_selection_index(
                    reference_ids,
                    current_value=cast(
                        str | None,
                        st.session_state.get(existing_reference_id_key),
                    ),
                    default_value=matched_reference_id or reference_ids[0],
                ),
                format_func=lambda reference_id: _reference_display_label(
                    next(
                        stored_reference
                        for stored_reference in stored_references
                        if stored_reference.reference_id == reference_id
                    )
                ),
                key=existing_reference_id_key,
            ),
        )
        stored_reference = get_stored_reference(
            st.session_state,
            reference_id=selected_reference_id,
        )
        if stored_reference is None:
            st.warning("The selected stored reference is no longer available.")
            return _ResolvedReferenceSelection(
                source_mode="existing",
                reference_name=None,
                reference_text=None,
                is_valid=False,
            )
        return _ResolvedReferenceSelection(
            source_mode="existing",
            reference_name=stored_reference.name,
            reference_text=stored_reference.reference_text,
            is_valid=True,
        )

    if not stored_references:
        st.caption("No stored references are available in this session yet.")

    input_mode_options: tuple[NewReferenceInputMode, ...] = ("upload", "paste")
    default_input_mode: NewReferenceInputMode = (
        "paste" if _normalized_optional_text(default_reference_text) else "upload"
    )
    input_mode = cast(
        NewReferenceInputMode,
        st.selectbox(
            "New reference input",
            options=list(input_mode_options),
            index=_selection_index(
                input_mode_options,
                current_value=cast(
                    str | None,
                    st.session_state.get(reference_input_mode_key),
                ),
                default_value=default_input_mode,
            ),
            format_func=lambda value: (
                "Upload file" if value == "upload" else "Paste text"
            ),
            key=reference_input_mode_key,
        ),
    )
    reference_name_value = st.text_input(
        "Reference name (optional)",
        value=default_reference_name or "",
        key=reference_name_key,
    )

    upload_col, paste_col = st.columns(2)
    with upload_col:
        uploaded_reference = st.file_uploader(
            "Upload FASTA/text reference",
            type=["fa", "fasta", "txt"],
            accept_multiple_files=False,
            key=reference_upload_key,
            disabled=input_mode != "upload",
        )
    with paste_col:
        pasted_reference_text = st.text_area(
            "Paste reference sequence",
            value=default_reference_text or "",
            key=reference_text_key,
            height=180,
            disabled=input_mode != "paste",
            placeholder="Paste a FASTA record or plain sequence here.",
        )

    raw_reference_text = (
        _uploaded_reference_text(uploaded_reference)
        if input_mode == "upload"
        else pasted_reference_text
    )
    if not raw_reference_text.strip():
        return _ResolvedReferenceSelection(
            source_mode="new",
            reference_name=_normalized_optional_text(reference_name_value),
            reference_text=None,
            is_valid=False,
        )

    try:
        resolved_reference_name, canonical_reference_text = _canonical_reference_input(
            reference_name=reference_name_value,
            reference_text=raw_reference_text,
            fallback_name=(
                _uploaded_reference_name(uploaded_reference)
                if input_mode == "upload"
                else None
            ),
        )
    except ValueError as exc:
        st.warning(str(exc))
        return _ResolvedReferenceSelection(
            source_mode="new",
            reference_name=_normalized_optional_text(reference_name_value),
            reference_text=None,
            is_valid=False,
        )

    st.caption(f"Resolved reference: {resolved_reference_name}")
    return _ResolvedReferenceSelection(
        source_mode="new",
        reference_name=resolved_reference_name,
        reference_text=canonical_reference_text,
        is_valid=True,
    )


def _render_pairwise_result(
    *,
    computed_alignment: ComputedAlignment,
    prepared_batch,
    theme_type: str,
) -> None:
    computed_assembly = computed_alignment.assembly
    if computed_assembly is None or computed_assembly.result is None:
        st.error(computed_alignment.status_reason or "Assembly result is unavailable.")
        return

    if not isinstance(computed_assembly.result, AssemblyResult):
        st.info("The new Alignments page currently renders only pairwise results.")
        return

    result = computed_assembly.result
    metric_col_1, metric_col_2, metric_col_3, metric_col_4, metric_col_5 = st.columns(5)
    with metric_col_1:
        st.metric(
            "Right orientation",
            result.chosen_right_orientation.replace("_", "-"),
        )
    with metric_col_2:
        st.metric(
            "Score",
            "NA" if result.score == float("-inf") else f"{result.score:.1f}",
        )
    with metric_col_3:
        st.metric("Overlap length", result.overlap_length)
    with metric_col_4:
        st.metric("Identity", f"{result.percent_identity:.1f}%")
    with metric_col_5:
        st.metric("Conflict columns", result.conflict_count)

    assembly_trace_view = build_selected_assembly_trace_view(
        prepared_batch,
        computed_assembly,
    )
    if assembly_trace_view is not None:
        st.subheader("Aligned electropherogram view")
        aligned_trace_figure = build_assembly_trace_figure(
            assembly_trace_view,
            theme_type=theme_type,
        )
        aligned_trace_figure.update_xaxes(
            range=_centered_alignment_x_range(
                alignment_length=assembly_trace_view.alignment_length,
                cell_width=assembly_trace_view.cell_width,
                center_column_index=None,
            )
        )
        aligned_trace_figure.update_layout(
            height=max(
                _ALIGNED_TRACE_MIN_HEIGHT_PX,
                _ALIGNED_TRACE_BASE_HEIGHT_PX
                + (len(assembly_trace_view.rows) * _ALIGNED_TRACE_ROW_HEIGHT_PX),
            ),
            margin={"l": 24, "r": 24, "t": 24, "b": 56},
        )
        st.plotly_chart(
            aligned_trace_figure,
            width="stretch",
            config={"scrollZoom": False},
        )

    st.subheader("Gapped alignment")
    st.code(format_assembly_block(result), wrap_lines=False)

    conflict_rows = assembly_conflicts_to_rows(result)
    if conflict_rows:
        st.subheader("Conflict columns")
        st.dataframe(conflict_rows, hide_index=True, width="stretch")
    else:
        st.success("No conflict columns were detected in this pairwise result.")

    if (
        computed_assembly.consensus_record is not None
        and computed_assembly.consensus_record.sequence
    ):
        consensus_record = computed_assembly.consensus_record
        consensus_fasta = to_fasta(consensus_record, line_width=None)
        alignment_fasta = format_assembly_alignment_fasta(
            result,
            consensus_name=consensus_record.name,
        )
        st.subheader("Consensus FASTA")
        st.code(consensus_fasta, wrap_lines=True)
        consensus_col, alignment_col = st.columns(2)
        filename_stem = _safe_filename_stem(consensus_record.name)
        with consensus_col:
            st.download_button(
                "Download consensus FASTA",
                data=consensus_fasta,
                file_name=f"{filename_stem}.fasta",
                mime="text/plain",
            )
        with alignment_col:
            st.download_button(
                "Download alignment FASTA",
                data=alignment_fasta,
                file_name=f"{filename_stem}.msa.fasta",
                mime="text/plain",
            )


def _render_reference_single_result(
    *,
    computed_alignment: ComputedAlignment,
    theme_type: str,
) -> None:
    reference_alignment = computed_alignment.reference_alignment
    if reference_alignment is None:
        st.error(
            computed_alignment.status_reason
            or "Reference-alignment result is unavailable."
        )
        return

    alignment_result = reference_alignment.alignment_result
    filename_stem = _safe_filename_stem(computed_alignment.definition.name)
    alignment_fasta = format_reference_alignment_fasta(alignment_result)
    event_rows = list(reference_alignment.event_rows)

    debug = False
    selected_event_row: dict[str, object] | None = None
    if debug and event_rows:
        selected_event_index = st.session_state.get(_SELECTED_EVENT_WIDGET_KEY, 0)
        if isinstance(selected_event_index, int) and 0 <= selected_event_index < len(
            event_rows
        ):
            selected_event_row = event_rows[selected_event_index]

    trace_view = reference_alignment.trace_view
    if trace_view is not None:
        selected_column_index = (
            None if selected_event_row is None else selected_event_row.get("column")
        )
        resolved_selected_column_index = (
            int(selected_column_index)
            if isinstance(selected_column_index, (int, float))
            else None
        )
        st.subheader("Aligned electropherogram view")
        aligned_trace_figure = build_reference_alignment_trace_figure(
            trace_view,
            theme_type=theme_type,
            selected_column_index=resolved_selected_column_index,
        )
        if resolved_selected_column_index is not None:
            aligned_trace_figure.update_xaxes(
                range=_centered_alignment_x_range(
                    alignment_length=trace_view.alignment_length,
                    cell_width=trace_view.cell_width,
                    center_column_index=resolved_selected_column_index,
                )
            )
        aligned_trace_figure.update_layout(
            height=max(
                _ALIGNED_TRACE_MIN_HEIGHT_PX,
                _ALIGNED_TRACE_BASE_HEIGHT_PX
                + (len(trace_view.rows) * _ALIGNED_TRACE_ROW_HEIGHT_PX),
            ),
            margin={"l": 24, "r": 24, "t": 24, "b": 56},
        )
        st.plotly_chart(
            aligned_trace_figure,
            width="stretch",
            config={"scrollZoom": False},
        )

    st.download_button(
        "Download gapped alignment FASTA",
        data=alignment_fasta,
        file_name=f"{filename_stem}.msa.fasta",
        mime="text/plain",
    )

    if debug:
        (
            metric_col_1,
            metric_col_2,
            metric_col_3,
            metric_col_4,
            metric_col_5,
            metric_col_6,
        ) = st.columns(6)
        with metric_col_1:
            st.metric("Strand", alignment_result.strand.replace("_", "-"))
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

        st.subheader("Alignment events")
        if event_rows:
            st.dataframe(event_rows, hide_index=True, width="stretch")
            st.selectbox(
                "Center chromatogram on event",
                options=list(range(len(event_rows))),
                format_func=lambda index: _event_option_label(event_rows[index]),
                key=_SELECTED_EVENT_WIDGET_KEY,
            )
        else:
            st.success(
                "No mismatch or indel events were detected in the current alignment."
            )


def _render_reference_multi_result(
    *,
    computed_alignment: ComputedAlignment,
    theme_type: str,
) -> None:
    reference_multi_alignment = computed_alignment.reference_multi_alignment
    if reference_multi_alignment is None:
        st.error(
            computed_alignment.status_reason
            or "Shared-reference alignment result is unavailable."
        )
        return

    result = reference_multi_alignment.result
    metric_col_1, metric_col_2, metric_col_3, metric_col_4, metric_col_5 = st.columns(5)
    with metric_col_1:
        st.metric("Included reads", result.included_member_count)
    with metric_col_2:
        st.metric("Excluded reads", result.excluded_member_count)
    with metric_col_3:
        st.metric("Reference length", len(result.reference_sequence))
    with metric_col_4:
        st.metric("Consensus length", len(result.consensus_sequence))
    with metric_col_5:
        st.metric("Ambiguous columns", result.ambiguous_column_count)

    trace_view = reference_multi_alignment.trace_view
    if trace_view is not None:
        st.subheader("Aligned electropherogram view")
        aligned_trace_figure = build_reference_multi_alignment_trace_figure(
            trace_view,
            theme_type=theme_type,
            selected_column_index=None,
        )
        aligned_trace_figure.update_layout(
            height=max(
                _ALIGNED_TRACE_MIN_HEIGHT_PX,
                _ALIGNED_TRACE_BASE_HEIGHT_PX
                + (len(trace_view.rows) * _ALIGNED_TRACE_ROW_HEIGHT_PX),
            ),
            margin={"l": 24, "r": 24, "t": 24, "b": 56},
        )
        st.plotly_chart(
            aligned_trace_figure,
            width="stretch",
            config={"scrollZoom": False},
        )
    else:
        st.warning(
            "Aligned electropherogram view is unavailable for this shared-reference alignment."
        )

    member_rows = list(reference_multi_alignment.member_rows)
    if member_rows:
        st.subheader("Read placement summary")
        st.dataframe(member_rows, hide_index=True, width="stretch")

    column_rows = list(reference_multi_alignment.column_rows)
    st.subheader("Support / variant summary")
    if column_rows:
        st.dataframe(column_rows, hide_index=True, width="stretch")
    else:
        st.success("No non-reference support columns remain after the current filters.")

    if result.consensus_sequence:
        consensus_record = consensus_record_from_reference_multi_result(
            result,
            name=computed_alignment.definition.name,
        )
        consensus_fasta = to_fasta(consensus_record, line_width=None)
        alignment_fasta = format_reference_multi_alignment_fasta(
            result,
            consensus_name=consensus_record.name,
        )
        st.subheader("Consensus FASTA")
        st.code(consensus_fasta, wrap_lines=True)
        consensus_col, alignment_col = st.columns(2)
        filename_stem = _safe_filename_stem(consensus_record.name)
        with consensus_col:
            st.download_button(
                "Download consensus FASTA",
                data=consensus_fasta,
                file_name=f"{filename_stem}.fasta",
                mime="text/plain",
            )
        with alignment_col:
            st.download_button(
                "Download alignment FASTA",
                data=alignment_fasta,
                file_name=f"{filename_stem}.msa.fasta",
                mime="text/plain",
            )

        debug = False
        if debug:
            st.subheader("Gapped alignment")
            st.code(
                format_reference_multi_alignment_block(result),
                wrap_lines=False,
            )


@st.dialog(
    "Create alignment",
    width="large",
    on_dismiss="rerun",
)
def _create_alignment_dialog(
    *,
    record_names: tuple[str, ...],
    record_annotations_by_name: dict[str, str],
    existing_names: tuple[str, ...],
) -> None:
    stored_references = list_stored_references(st.session_state)
    default_reads = tuple(record_names[:2])
    default_engine_kind: AlignmentEngineKind = "pairwise"
    alignment_name = st.text_input(
        "Alignment name",
        value=_default_dialog_name(
            source_filenames=default_reads,
            engine_kind=default_engine_kind,
            reference_name=None,
            existing_names=existing_names,
        ),
        key=_NEW_NAME_WIDGET_KEY,
    )
    engine_kind = cast(
        AlignmentEngineKind,
        st.selectbox(
            "Engine",
            options=["pairwise", "reference_single", "reference_multi"],
            format_func=_engine_label,
            key=_NEW_ENGINE_WIDGET_KEY,
        ),
    )
    selected_reads = cast(
        list[str],
        st.multiselect(
            "Reads",
            options=record_names,
            default=list(default_reads),
            format_func=lambda filename: record_annotations_by_name[filename],
            key=_NEW_READS_WIDGET_KEY,
        ),
    )
    reference_selection = _render_reference_selection(
        engine_kind=engine_kind,
        stored_references=stored_references,
        reference_source_key=_NEW_REFERENCE_SOURCE_WIDGET_KEY,
        existing_reference_id_key=_NEW_EXISTING_REFERENCE_ID_WIDGET_KEY,
        reference_input_mode_key=_NEW_REFERENCE_INPUT_MODE_WIDGET_KEY,
        reference_name_key=_NEW_REFERENCE_NAME_WIDGET_KEY,
        reference_text_key=_NEW_REFERENCE_TEXT_WIDGET_KEY,
        reference_upload_key=_NEW_REFERENCE_UPLOAD_WIDGET_KEY,
    )
    strand_policy = cast(
        StrandPolicy,
        st.selectbox(
            "Strand policy",
            options=["auto", "forward", "reverse_complement"],
            format_func=lambda value: value.replace("_", "-"),
            key=_NEW_STRAND_POLICY_WIDGET_KEY,
            disabled=not _reference_required(engine_kind),
        ),
    )

    min_overlap_length = int(
        st.number_input(
            "Minimum overlap length",
            min_value=1,
            step=1,
            value=25,
            key=_NEW_MIN_OVERLAP_WIDGET_KEY,
            disabled=engine_kind == "reference_single",
        )
    )
    min_percent_identity = float(
        st.number_input(
            "Minimum percent identity",
            min_value=0.0,
            max_value=100.0,
            step=1.0,
            value=70.0,
            key=_NEW_MIN_IDENTITY_WIDGET_KEY,
            disabled=engine_kind == "reference_single",
        )
    )
    quality_margin = int(
        st.number_input(
            "Quality margin for mismatch resolution",
            min_value=1,
            step=1,
            value=3,
            key=_NEW_QUALITY_MARGIN_WIDGET_KEY,
            disabled=engine_kind == "reference_single",
        )
    )

    selection_valid = _alignment_member_count_valid(
        engine_kind,
        selected_read_count=len(selected_reads),
    )
    reference_valid = (not _reference_required(engine_kind)) or (
        reference_selection.is_valid
        and reference_selection.reference_name is not None
        and reference_selection.reference_text is not None
    )
    if not selection_valid:
        if engine_kind == "pairwise":
            st.warning(
                "Pairwise alignment currently requires exactly 2 selected reads."
            )
        elif engine_kind == "reference_single":
            st.warning(
                "Single-reference alignment currently requires exactly 1 selected read."
            )
        else:
            st.warning(
                "Reference-guided multi-read alignment currently requires at least 2 reads."
            )
    elif _reference_required(engine_kind) and not reference_valid:
        st.warning("Reference-guided modes require one valid reference.")
    elif engine_kind == "reference_multi":
        st.caption(
            "Shared-reference mode independently places each selected read onto "
            "the reference coordinate system, then derives support and consensus "
            "across the shared grid."
        )

    save_col, cancel_col = st.columns(2)
    with save_col:
        save_clicked = st.button(
            "Create alignment",
            type="primary",
            use_container_width=True,
            disabled=not (selection_valid and reference_valid),
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
        )

    if cancel_clicked:
        st.rerun()

    if save_clicked:
        reference_name = reference_selection.reference_name
        reference_text = reference_selection.reference_text
        if _reference_required(engine_kind):
            if reference_name is None or reference_text is None:
                st.warning("Reference-guided modes require one valid reference.")
                return
            if reference_selection.source_mode == "new":
                stored_reference = store_reference(
                    st.session_state,
                    name=reference_name,
                    reference_text=reference_text,
                )
                reference_name = stored_reference.name
                reference_text = stored_reference.reference_text

        create_alignment_definition(
            st.session_state,
            name=alignment_name,
            source_filenames=selected_reads,
            engine_kind=engine_kind,
            assembly_config=AssemblyConfig(
                min_overlap_length=min_overlap_length,
                min_percent_identity=min_percent_identity,
                quality_margin=quality_margin,
            ),
            reference_name=reference_name,
            reference_text=reference_text,
            strand_policy=strand_policy,
        )
        st.rerun()


@st.dialog(
    "Edit alignment",
    width="large",
    on_dismiss="rerun",
)
def _edit_alignment_dialog(
    definition: AlignmentDefinition,
    *,
    record_names: tuple[str, ...],
    record_annotations_by_name: dict[str, str],
    existing_names: tuple[str, ...],
) -> None:
    stored_references = list_stored_references(st.session_state)
    alignment_name = st.text_input(
        "Alignment name",
        value=definition.name,
        key=f"{_EDIT_NAME_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
    )
    engine_kind = cast(
        AlignmentEngineKind,
        st.selectbox(
            "Engine",
            options=["pairwise", "reference_single", "reference_multi"],
            index=(
                0
                if definition.engine_kind == "pairwise"
                else (1 if definition.engine_kind == "reference_single" else 2)
            ),
            format_func=_engine_label,
            key=f"{_EDIT_ENGINE_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
        ),
    )
    selected_reads = cast(
        list[str],
        st.multiselect(
            "Reads",
            options=record_names,
            default=list(definition.source_filenames),
            format_func=lambda filename: record_annotations_by_name[filename],
            key=f"{_EDIT_READS_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
        ),
    )
    reference_selection = _render_reference_selection(
        engine_kind=engine_kind,
        stored_references=stored_references,
        reference_source_key=(
            f"{_EDIT_REFERENCE_SOURCE_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        existing_reference_id_key=(
            f"{_EDIT_EXISTING_REFERENCE_ID_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        reference_input_mode_key=(
            f"{_EDIT_REFERENCE_INPUT_MODE_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        reference_name_key=(
            f"{_EDIT_REFERENCE_NAME_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        reference_text_key=(
            f"{_EDIT_REFERENCE_TEXT_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        reference_upload_key=(
            f"{_EDIT_REFERENCE_UPLOAD_WIDGET_KEY_PREFIX}.{definition.alignment_id}"
        ),
        default_reference_name=definition.reference_name,
        default_reference_text=definition.reference_text,
    )
    strand_policy = cast(
        StrandPolicy,
        st.selectbox(
            "Strand policy",
            options=["auto", "forward", "reverse_complement"],
            index=(
                0
                if definition.strand_policy == "auto"
                else (1 if definition.strand_policy == "forward" else 2)
            ),
            format_func=lambda value: value.replace("_", "-"),
            key=f"{_EDIT_STRAND_POLICY_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
            disabled=not _reference_required(engine_kind),
        ),
    )

    min_overlap_length = int(
        st.number_input(
            "Minimum overlap length",
            min_value=1,
            step=1,
            value=definition.assembly_config.min_overlap_length,
            key=f"{_EDIT_MIN_OVERLAP_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
            disabled=engine_kind == "reference_single",
        )
    )
    min_percent_identity = float(
        st.number_input(
            "Minimum percent identity",
            min_value=0.0,
            max_value=100.0,
            step=1.0,
            value=definition.assembly_config.min_percent_identity,
            key=f"{_EDIT_MIN_IDENTITY_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
            disabled=engine_kind == "reference_single",
        )
    )
    quality_margin = int(
        st.number_input(
            "Quality margin for mismatch resolution",
            min_value=1,
            step=1,
            value=definition.assembly_config.quality_margin,
            key=f"{_EDIT_QUALITY_MARGIN_WIDGET_KEY_PREFIX}.{definition.alignment_id}",
            disabled=engine_kind == "reference_single",
        )
    )

    selection_valid = _alignment_member_count_valid(
        engine_kind,
        selected_read_count=len(selected_reads),
    )
    reference_valid = (not _reference_required(engine_kind)) or (
        reference_selection.is_valid
        and reference_selection.reference_name is not None
        and reference_selection.reference_text is not None
    )
    if not selection_valid:
        if engine_kind == "pairwise":
            st.warning(
                "Pairwise alignment currently requires exactly 2 selected reads."
            )
        elif engine_kind == "reference_single":
            st.warning(
                "Single-reference alignment currently requires exactly 1 selected read."
            )
        else:
            st.warning(
                "Reference-guided multi-read alignment currently requires at least 2 reads."
            )
    elif _reference_required(engine_kind) and not reference_valid:
        st.warning("Reference-guided modes require one valid reference.")
    elif engine_kind == "reference_multi":
        st.caption(
            "Shared-reference mode independently places each selected read onto "
            "the reference coordinate system, then derives support and consensus "
            "across the shared grid."
        )

    save_col, cancel_col = st.columns(2)
    with save_col:
        save_clicked = st.button(
            "Save changes",
            type="primary",
            use_container_width=True,
            disabled=not (selection_valid and reference_valid),
            key=f"alignments.dialog.edit.save.{definition.alignment_id}",
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
            key=f"alignments.dialog.edit.cancel.{definition.alignment_id}",
        )

    if cancel_clicked:
        st.rerun()

    if save_clicked:
        reference_name = reference_selection.reference_name
        reference_text = reference_selection.reference_text
        if _reference_required(engine_kind):
            if reference_name is None or reference_text is None:
                st.warning("Reference-guided modes require one valid reference.")
                return
            if reference_selection.source_mode == "new":
                stored_reference = store_reference(
                    st.session_state,
                    name=reference_name,
                    reference_text=reference_text,
                )
                reference_name = stored_reference.name
                reference_text = stored_reference.reference_text

        update_alignment_definition(
            st.session_state,
            alignment_id=definition.alignment_id,
            name=alignment_name,
            source_filenames=selected_reads,
            engine_kind=engine_kind,
            assembly_config=AssemblyConfig(
                min_overlap_length=min_overlap_length,
                min_percent_identity=min_percent_identity,
                quality_margin=quality_margin,
                match_score=definition.assembly_config.match_score,
                mismatch_score=definition.assembly_config.mismatch_score,
                open_internal_gap_score=definition.assembly_config.open_internal_gap_score,
                extend_internal_gap_score=definition.assembly_config.extend_internal_gap_score,
            ),
            reference_name=reference_name,
            reference_text=reference_text,
            strand_policy=strand_policy,
        )
        st.rerun()


@st.dialog(
    "Delete alignment",
    width="small",
    on_dismiss="rerun",
)
def _delete_alignment_dialog(definition: AlignmentDefinition) -> None:
    st.warning(
        f"Delete the saved alignment '{definition.name}'? "
        "This removes only the saved definition, not the source samples."
    )

    confirm_col, cancel_col = st.columns(2)
    with confirm_col:
        confirm_clicked = st.button(
            "Delete",
            type="primary",
            use_container_width=True,
            key=f"alignments.dialog.delete.confirm.{definition.alignment_id}",
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
            key=f"alignments.dialog.delete.cancel.{definition.alignment_id}",
        )

    if cancel_clicked:
        st.rerun()

    if confirm_clicked:
        delete_alignment_definition(
            st.session_state,
            alignment_id=definition.alignment_id,
        )
        st.rerun()


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
record_names = tuple(parsed_records.keys())
record_annotations = build_record_annotations(
    parsed_records.keys(),
    trim_state.trim_configs_by_record,
)
alignment_state = sync_alignment_session_state(
    st.session_state,
    batch_signature=parsed_batch.signature,
    parsed_record_names=record_names,
)
alignment_definitions_by_id = alignment_state.alignments_by_id

if not alignment_definitions_by_id:
    st.info(
        "No saved alignments are defined for the current batch yet. "
        "Create one to keep a reusable read/reference specification."
    )
    if st.button("Create alignment", type="primary"):
        _create_alignment_dialog(
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=(),
        )
    st.stop()

prepared_batch = prepare_batch_for_trim_state(
    parsed_batch,
    get_batch_trim_state(st.session_state),
)
computed_alignments = compute_saved_alignments(
    prepared_batch,
    tuple(alignment_definitions_by_id.values()),
)

selected_alignment_id = alignment_state.selected_alignment_id
select_options: list[str] = list(alignment_definitions_by_id)
select_index = (
    select_options.index(selected_alignment_id)
    if selected_alignment_id in alignment_definitions_by_id
    else 0
)

select_col, new_col, edit_col, delete_col = st.columns([5, 1, 1, 1])
with select_col:
    selected_alignment_id = cast(
        str,
        st.selectbox(
            "Alignment",
            options=select_options,
            index=select_index,
            format_func=lambda alignment_id: _alignment_select_label(
                alignment_definitions_by_id,
                alignment_id,
            ),
            key=_ALIGNMENT_SELECT_WIDGET_KEY,
        ),
    )
set_selected_alignment_id(st.session_state, selected_alignment_id)

selected_definition = alignment_definitions_by_id[selected_alignment_id]

with new_col:
    if st.button("New", use_container_width=True):
        _create_alignment_dialog(
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=tuple(
                definition.name for definition in alignment_definitions_by_id.values()
            ),
        )

with edit_col:
    if st.button("Edit", use_container_width=True):
        _edit_alignment_dialog(
            selected_definition,
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=tuple(
                definition.name
                for definition in alignment_definitions_by_id.values()
                if definition.alignment_id != selected_definition.alignment_id
            ),
        )

with delete_col:
    if st.button("Delete", use_container_width=True):
        _delete_alignment_dialog(selected_definition)

selected_computed_alignment = computed_alignments[selected_definition.alignment_id]

status_message = selected_computed_alignment.status_reason
if selected_computed_alignment.status == "ok":
    st.success(
        "Saved alignment resolved successfully against the current trimmed batch."
    )
elif selected_computed_alignment.status == "rejected":
    st.warning(status_message or "Saved alignment did not pass the current filters.")
elif selected_computed_alignment.status == "not_implemented":
    st.info(status_message or "Saved alignment engine is not implemented yet.")
else:
    st.error(status_message or "Saved alignment definition could not be resolved.")

debug = False
if debug:
    summary: dict[str, object] = {
        "engine_kind": selected_definition.engine_kind,
        "members": list(selected_definition.source_filenames),
    }
    if _reference_required(selected_definition.engine_kind):
        summary["reference_name"] = selected_definition.reference_name
        summary["strand_policy"] = selected_definition.strand_policy
    st.write(summary)

theme_type = str(getattr(getattr(st.context, "theme", None), "type", "light"))
if selected_definition.engine_kind == "pairwise":
    _render_pairwise_result(
        computed_alignment=selected_computed_alignment,
        prepared_batch=prepared_batch,
        theme_type=theme_type,
    )
elif selected_definition.engine_kind == "reference_single":
    _render_reference_single_result(
        computed_alignment=selected_computed_alignment,
        theme_type=theme_type,
    )
elif selected_definition.engine_kind == "reference_multi":
    _render_reference_multi_result(
        computed_alignment=selected_computed_alignment,
        theme_type=theme_type,
    )
else:
    st.error("Unsupported alignment engine.")
