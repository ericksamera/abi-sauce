from __future__ import annotations

from typing import cast
from urllib.parse import quote

import streamlit as st

from abi_sauce.assembly import (
    AssemblyConfig,
    assembly_conflicts_to_rows,
    format_assembly_block,
)
from abi_sauce.assembly_state import (
    AssemblyDefinition,
    create_assembly_definition,
    delete_assembly_definition,
    set_export_selected_assembly_ids,
    set_selected_assembly_id,
    suggest_assembly_name,
    sync_assembly_session_state,
    update_assembly_definition,
)
from abi_sauce.chromatogram import ChromatogramView, build_chromatogram_view
from abi_sauce.chromatogram_figure import build_chromatogram_figure
from abi_sauce.export import to_fasta
from abi_sauce.exceptions import ExportError
from abi_sauce.services.assembly import (
    ComputedAssembly,
    compute_saved_assemblies,
    prepare_assembly_download,
)
from abi_sauce.services.batch import apply_trim_configs
from abi_sauce.trim_state import build_record_annotations, resolve_batch_trim_inputs
from abi_sauce.upload_state import get_active_parsed_batch
from abi_sauce.viewer_state import get_batch_trim_state, sync_viewer_session_state

_ASSEMBLY_SELECT_WIDGET_KEY = "assembly.selected_assembly_widget"
_SELECTED_CONFLICT_WIDGET_KEY = "assembly.selected_conflict"
_EXPORT_SELECTED_WIDGET_KEY = "assembly.export_selected_widget"
_EXPORT_CONCATENATE_WIDGET_KEY = "assembly.export_concatenate"
_EXPORT_INCLUDE_REJECTED_WIDGET_KEY = "assembly.export_include_rejected"
_EXPORT_INCLUDE_MANIFEST_WIDGET_KEY = "assembly.export_include_manifest"
_EXPORT_FILENAME_STEM_WIDGET_KEY = "assembly.export_filename_stem"
_ASSEMBLY_SIDEBAR_EXPORT_POPOVER_KEY = "assembly.sidebar.export_popover"
_ASSEMBLY_SIDEBAR_EXPORT_CONCATENATE_WIDGET_KEY = "assembly.sidebar.export_concatenate"
_ASSEMBLY_SIDEBAR_EXPORT_INCLUDE_MANIFEST_WIDGET_KEY = (
    "assembly.sidebar.export_include_manifest"
)
_ASSEMBLY_SIDEBAR_EXPORT_FILENAME_STEM_WIDGET_KEY = (
    "assembly.sidebar.export_filename_stem"
)
_AUTO_PROMPTED_SIGNATURE_KEY = "abi_sauce.assembly.auto_prompted_signature"

_NEW_NAME_WIDGET_KEY = "assembly.dialog.new.name"
_NEW_READS_WIDGET_KEY = "assembly.dialog.new.reads"
_NEW_MIN_OVERLAP_WIDGET_KEY = "assembly.dialog.new.min_overlap"
_NEW_MIN_IDENTITY_WIDGET_KEY = "assembly.dialog.new.min_identity"
_NEW_QUALITY_MARGIN_WIDGET_KEY = "assembly.dialog.new.quality_margin"

_EDIT_NAME_WIDGET_KEY_PREFIX = "assembly.dialog.edit.name"
_EDIT_READS_WIDGET_KEY_PREFIX = "assembly.dialog.edit.reads"
_EDIT_MIN_OVERLAP_WIDGET_KEY_PREFIX = "assembly.dialog.edit.min_overlap"
_EDIT_MIN_IDENTITY_WIDGET_KEY_PREFIX = "assembly.dialog.edit.min_identity"
_EDIT_QUALITY_MARGIN_WIDGET_KEY_PREFIX = "assembly.dialog.edit.quality_margin"

st.set_page_config(page_title="Assembly", layout="wide")
st.title("Assembly")
st.caption(
    "Manage saved assembly definitions, inspect their recomputed results "
    "against the current trimmed batch, and export consensus sequences."
)


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


def _conflict_option_label(conflict_row: dict[str, object]) -> str:
    return (
        f"col {conflict_row['column']} | "
        f"{conflict_row['left_base']}>{conflict_row['right_base']} | "
        f"{conflict_row['resolution']}"
    )


def _assembly_select_label(
    assembly_definitions_by_id: dict[str, AssemblyDefinition],
    assembly_id: str,
) -> str:
    definition = assembly_definitions_by_id[assembly_id]
    return definition.name


def _default_dialog_name(
    *,
    source_filenames: tuple[str, ...],
    existing_names: tuple[str, ...],
    current_name: str | None = None,
) -> str:
    if current_name is not None and current_name.strip():
        return current_name
    return suggest_assembly_name(
        source_filenames,
        existing_names=existing_names,
    )


def _assembly_summary_rows(
    computed_assemblies: dict[str, ComputedAssembly],
    *,
    export_selected_ids: frozenset[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for assembly_id, computed_assembly in computed_assemblies.items():
        result = computed_assembly.result
        rows.append(
            {
                "assembly": computed_assembly.definition.name,
                "members": ", ".join(computed_assembly.definition.source_filenames),
                "selected_for_export": assembly_id in export_selected_ids,
                "status": computed_assembly.status,
                "status_reason": computed_assembly.status_reason,
                "right_orientation": (
                    None if result is None else result.chosen_right_orientation
                ),
                "overlap_length": None if result is None else result.overlap_length,
                "percent_identity": None if result is None else result.percent_identity,
                "conflict_count": None if result is None else result.conflict_count,
            }
        )
    return rows


def _effective_export_selected_ids(
    assembly_definitions_by_id: dict[str, AssemblyDefinition],
    export_selected_ids: frozenset[str],
) -> frozenset[str]:
    return frozenset(assembly_definitions_by_id)


def _build_assembly_blast_url(
    computed_assemblies: dict[str, ComputedAssembly],
    *,
    selected_ids: frozenset[str],
) -> str | None:
    download_artifact = prepare_assembly_download(
        computed_assemblies,
        selected_ids=selected_ids,
        concatenate_batch=True,
        filename_stem="abi-sauce-assemblies",
        require_accepted=True,
        include_manifest=False,
        fasta_line_width=None,
    )
    if not download_artifact.is_downloadable or not isinstance(
        download_artifact.data, str
    ):
        return None

    query_text = download_artifact.data.strip()
    if not query_text:
        return None

    return (
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        "?PROGRAM=blastn"
        "&PAGE_TYPE=BlastSearch"
        "&LINK_LOC=blasthome"
        f"&QUERY={quote(query_text, safe='')}"
    )


def _render_assembly_sidebar_controls(
    *,
    assembly_definitions_by_id: dict[str, AssemblyDefinition],
    computed_assemblies: dict[str, ComputedAssembly],
    export_selected_ids: frozenset[str],
) -> frozenset[str]:
    effective_export_selected_ids = _effective_export_selected_ids(
        assembly_definitions_by_id,
        export_selected_ids,
    )
    if effective_export_selected_ids != export_selected_ids:
        set_export_selected_assembly_ids(
            st.session_state,
            effective_export_selected_ids,
        )

    with st.sidebar:
        with st.popover(
            "Export Assemblies",
            key=_ASSEMBLY_SIDEBAR_EXPORT_POPOVER_KEY,
            width="stretch",
            icon=":material/download:",
        ):
            st.caption("Exports all accepted saved assemblies from the current batch.")

            concatenate_batch = st.checkbox(
                "Concatenate entries into a single FASTA",
                value=True,
                key=_ASSEMBLY_SIDEBAR_EXPORT_CONCATENATE_WIDGET_KEY,
            )
            if not concatenate_batch:
                st.caption(
                    "Individual consensus records will be bundled into a ZIP archive."
                )

            include_manifest = st.checkbox(
                "Include manifest.json in ZIP export",
                value=True,
                key=_ASSEMBLY_SIDEBAR_EXPORT_INCLUDE_MANIFEST_WIDGET_KEY,
                disabled=concatenate_batch,
            )
            filename_stem = st.text_input(
                "Output filename stem",
                value="abi-sauce-assemblies",
                key=_ASSEMBLY_SIDEBAR_EXPORT_FILENAME_STEM_WIDGET_KEY,
            )

            try:
                download_artifact = prepare_assembly_download(
                    computed_assemblies,
                    selected_ids=effective_export_selected_ids,
                    concatenate_batch=bool(concatenate_batch),
                    filename_stem=filename_stem,
                    require_accepted=True,
                    include_manifest=bool(include_manifest),
                    fasta_line_width=None,
                )
            except ExportError as exc:
                st.warning(str(exc))
            else:
                st.caption(
                    f"{len(download_artifact.eligible_records)} exportable consensus "
                    f"record(s) · {len(download_artifact.ineligible_reasons)} excluded"
                )
                if download_artifact.ineligible_reasons:
                    with st.expander("Excluded from current export"):
                        for (
                            assembly_name,
                            reasons,
                        ) in download_artifact.ineligible_reasons:
                            st.warning(f"{assembly_name}: {', '.join(reasons)}")

                if not download_artifact.is_downloadable:
                    st.warning(
                        "No saved assemblies are eligible for this export selection."
                    )
                else:
                    st.download_button(
                        "Export Assemblies",
                        data=download_artifact.data,
                        file_name=download_artifact.filename,
                        mime=download_artifact.mime,
                        key="assembly.sidebar.export_button",
                    )

        st.divider()

        blast_url = _build_assembly_blast_url(
            computed_assemblies,
            selected_ids=effective_export_selected_ids,
        )
        if blast_url is None:
            st.caption(
                "No accepted assembly consensus sequences are available for BLAST."
            )
        elif len(blast_url) > 60_000:
            st.caption(
                "Assembly batch is too large for a reliable BLAST URL; export FASTA "
                "and paste/upload it in BLAST instead."
            )
        else:
            st.link_button(
                "BLAST Assemblies",
                blast_url,
                width="stretch",
                icon=":material/rocket_launch:",
            )

    return effective_export_selected_ids


@st.dialog(
    "Create assembly",
    width="small",
    on_dismiss="rerun",
)
def _create_assembly_dialog(
    *,
    record_names: tuple[str, ...],
    record_annotations_by_name: dict[str, str],
    existing_names: tuple[str, ...],
) -> None:
    st.caption(
        "Saved assemblies store only the read selection and thresholds. "
        "The actual assembly result is recomputed from the current trimmed batch "
        "on each rerun."
    )

    default_reads = tuple(record_names[:2])
    suggested_name = _default_dialog_name(
        source_filenames=default_reads,
        existing_names=existing_names,
    )
    assembly_name = st.text_input(
        "Assembly name",
        value=suggested_name,
        key=_NEW_NAME_WIDGET_KEY,
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
    min_overlap_length = int(
        st.number_input(
            "Minimum overlap length",
            min_value=1,
            step=1,
            value=25,
            key=_NEW_MIN_OVERLAP_WIDGET_KEY,
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
        )
    )
    quality_margin = int(
        st.number_input(
            "Quality margin for mismatch resolution",
            min_value=1,
            step=1,
            value=3,
            key=_NEW_QUALITY_MARGIN_WIDGET_KEY,
        )
    )

    pairwise_valid = len(selected_reads) == 2
    if not pairwise_valid:
        st.warning("Current pairwise engine requires exactly 2 selected reads.")

    save_col, cancel_col = st.columns(2)
    with save_col:
        save_clicked = st.button(
            "Create assembly",
            type="primary",
            use_container_width=True,
            disabled=not pairwise_valid,
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
        )

    if cancel_clicked:
        st.rerun()

    if save_clicked:
        create_assembly_definition(
            st.session_state,
            name=assembly_name,
            source_filenames=selected_reads,
            config=AssemblyConfig(
                min_overlap_length=min_overlap_length,
                min_percent_identity=min_percent_identity,
                quality_margin=quality_margin,
            ),
        )
        st.rerun()


@st.dialog(
    "Edit assembly",
    width="small",
    on_dismiss="rerun",
)
def _edit_assembly_dialog(
    definition: AssemblyDefinition,
    *,
    record_names: tuple[str, ...],
    record_annotations_by_name: dict[str, str],
    existing_names: tuple[str, ...],
) -> None:
    assembly_name = st.text_input(
        "Assembly name",
        value=definition.name,
        key=f"{_EDIT_NAME_WIDGET_KEY_PREFIX}.{definition.assembly_id}",
    )
    selected_reads = cast(
        list[str],
        st.multiselect(
            "Reads",
            options=record_names,
            default=list(definition.source_filenames),
            format_func=lambda filename: record_annotations_by_name[filename],
            key=f"{_EDIT_READS_WIDGET_KEY_PREFIX}.{definition.assembly_id}",
        ),
    )
    min_overlap_length = int(
        st.number_input(
            "Minimum overlap length",
            min_value=1,
            step=1,
            value=definition.config.min_overlap_length,
            key=f"{_EDIT_MIN_OVERLAP_WIDGET_KEY_PREFIX}.{definition.assembly_id}",
        )
    )
    min_percent_identity = float(
        st.number_input(
            "Minimum percent identity",
            min_value=0.0,
            max_value=100.0,
            step=1.0,
            value=definition.config.min_percent_identity,
            key=f"{_EDIT_MIN_IDENTITY_WIDGET_KEY_PREFIX}.{definition.assembly_id}",
        )
    )
    quality_margin = int(
        st.number_input(
            "Quality margin for mismatch resolution",
            min_value=1,
            step=1,
            value=definition.config.quality_margin,
            key=f"{_EDIT_QUALITY_MARGIN_WIDGET_KEY_PREFIX}.{definition.assembly_id}",
        )
    )

    pairwise_valid = len(selected_reads) == 2
    if not pairwise_valid:
        st.warning("Current pairwise engine requires exactly 2 selected reads.")

    save_col, cancel_col = st.columns(2)
    with save_col:
        save_clicked = st.button(
            "Save changes",
            type="primary",
            use_container_width=True,
            disabled=not pairwise_valid,
            key=f"assembly.dialog.edit.save.{definition.assembly_id}",
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
            key=f"assembly.dialog.edit.cancel.{definition.assembly_id}",
        )

    if cancel_clicked:
        st.rerun()

    if save_clicked:
        update_assembly_definition(
            st.session_state,
            assembly_id=definition.assembly_id,
            name=assembly_name,
            source_filenames=selected_reads,
            config=AssemblyConfig(
                min_overlap_length=min_overlap_length,
                min_percent_identity=min_percent_identity,
                quality_margin=quality_margin,
                match_score=definition.config.match_score,
                mismatch_score=definition.config.mismatch_score,
                open_internal_gap_score=definition.config.open_internal_gap_score,
                extend_internal_gap_score=definition.config.extend_internal_gap_score,
            ),
        )
        st.rerun()


@st.dialog(
    "Delete assembly",
    width="small",
    on_dismiss="rerun",
)
def _delete_assembly_dialog(definition: AssemblyDefinition) -> None:
    st.warning(
        f"Delete the saved assembly '{definition.name}'? "
        "This removes only the saved definition, not the source samples."
    )

    confirm_col, cancel_col = st.columns(2)
    with confirm_col:
        confirm_clicked = st.button(
            "Delete",
            type="primary",
            use_container_width=True,
            key=f"assembly.dialog.delete.confirm.{definition.assembly_id}",
        )
    with cancel_col:
        cancel_clicked = st.button(
            "Cancel",
            use_container_width=True,
            key=f"assembly.dialog.delete.cancel.{definition.assembly_id}",
        )

    if cancel_clicked:
        st.rerun()

    if confirm_clicked:
        delete_assembly_definition(
            st.session_state,
            assembly_id=definition.assembly_id,
        )
        st.rerun()


parsed_batch = get_active_parsed_batch(st.session_state)
if parsed_batch is None:
    st.info("Load one or more .ab1 files from the workspace sidebar to begin.")
    st.stop()

parsed_records = parsed_batch.parsed_records
if len(parsed_records) < 2:
    st.warning("Assembly needs at least two parsed records in the active batch.")
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
assembly_state = sync_assembly_session_state(
    st.session_state,
    batch_signature=parsed_batch.signature,
    parsed_record_names=record_names,
)
assembly_definitions_by_id = assembly_state.assemblies_by_id

if not assembly_definitions_by_id:
    auto_prompted_signature = st.session_state.get(_AUTO_PROMPTED_SIGNATURE_KEY)
    if auto_prompted_signature != parsed_batch.signature:
        st.session_state[_AUTO_PROMPTED_SIGNATURE_KEY] = parsed_batch.signature
        _create_assembly_dialog(
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=(),
        )

    st.info(
        "No saved assemblies are defined for the current batch yet. "
        "Create one to keep a reusable read pairing and exportable consensus."
    )
    if st.button("Create assembly", type="primary"):
        _create_assembly_dialog(
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=(),
        )
    st.stop()

resolved_trim_inputs = resolve_batch_trim_inputs(get_batch_trim_state(st.session_state))
prepared_batch = apply_trim_configs(
    parsed_batch,
    default_trim_config=resolved_trim_inputs.default_trim_config,
    trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
)
computed_assemblies = compute_saved_assemblies(
    prepared_batch,
    tuple(assembly_definitions_by_id.values()),
)
effective_export_selected_ids = _render_assembly_sidebar_controls(
    assembly_definitions_by_id=assembly_definitions_by_id,
    computed_assemblies=computed_assemblies,
    export_selected_ids=assembly_state.export_selected_assembly_ids,
)

selected_assembly_id = assembly_state.selected_assembly_id
select_options: list[str] = list(assembly_definitions_by_id)
select_index = (
    select_options.index(selected_assembly_id)
    if selected_assembly_id in assembly_definitions_by_id
    else 0
)

select_col, new_col, edit_col, delete_col = st.columns([5, 1, 1, 1])
with select_col:
    selected_assembly_id = cast(
        str,
        st.selectbox(
            "Assembly",
            options=select_options,
            index=select_index,
            format_func=lambda assembly_id: _assembly_select_label(
                assembly_definitions_by_id,
                assembly_id,
            ),
            key=_ASSEMBLY_SELECT_WIDGET_KEY,
        ),
    )
set_selected_assembly_id(st.session_state, selected_assembly_id)

selected_definition = assembly_definitions_by_id[selected_assembly_id]

with new_col:
    if st.button("New", use_container_width=True):
        _create_assembly_dialog(
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=tuple(
                definition.name for definition in assembly_definitions_by_id.values()
            ),
        )

with edit_col:
    if (
        st.button(
            "Edit",
            use_container_width=True,
            disabled=selected_definition is None,
        )
        and selected_definition is not None
    ):
        _edit_assembly_dialog(
            selected_definition,
            record_names=record_names,
            record_annotations_by_name=record_annotations.display_labels_by_record,
            existing_names=tuple(
                definition.name
                for definition in assembly_definitions_by_id.values()
                if definition.assembly_id != selected_definition.assembly_id
            ),
        )

with delete_col:
    if (
        st.button(
            "Delete",
            use_container_width=True,
            disabled=selected_definition is None,
        )
        and selected_definition is not None
    ):
        _delete_assembly_dialog(selected_definition)

st.subheader("Saved assemblies")
st.dataframe(
    _assembly_summary_rows(
        computed_assemblies,
        export_selected_ids=effective_export_selected_ids,
    ),
    hide_index=True,
    width="stretch",
)


if selected_definition is None:
    st.info("Select a saved assembly to inspect its current result.")
    st.stop()

selected_computed_assembly = computed_assemblies[selected_definition.assembly_id]
selected_result = selected_computed_assembly.result

member_rows = [
    {
        "member": f"read {index}",
        "filename": source_filename,
        "display_name": prepared_batch.parsed_records[source_filename].name,
        "display_orientation": prepared_batch.parsed_records[
            source_filename
        ].orientation,
        "trimmed_length": prepared_batch.trim_results[source_filename].trimmed_length,
        "has_qualities": prepared_batch.trim_results[source_filename].record.qualities
        is not None,
        "has_trace_data": prepared_batch.parsed_records[source_filename].trace_data
        is not None,
    }
    for index, source_filename in enumerate(
        selected_definition.source_filenames,
        start=1,
    )
]
st.subheader("Assembly members")
st.dataframe(member_rows, hide_index=True, width="stretch")

if selected_computed_assembly.status == "ok":
    st.success("Assembly candidate passed the current overlap and identity thresholds.")
elif selected_computed_assembly.status == "rejected":
    st.warning(
        selected_computed_assembly.status_reason
        or "Assembly candidate did not pass the current filters."
    )
else:
    st.error(
        selected_computed_assembly.status_reason
        or "Saved assembly definition could not be resolved."
    )

if selected_result is None:
    st.stop()

metric_col_1, metric_col_2, metric_col_3, metric_col_4, metric_col_5 = st.columns(5)
with metric_col_1:
    st.metric("Right orientation", selected_result.chosen_right_orientation)
with metric_col_2:
    st.metric(
        "Score",
        (
            "NA"
            if selected_result.score == float("-inf")
            else f"{selected_result.score:.1f}"
        ),
    )
with metric_col_3:
    st.metric("Overlap length", selected_result.overlap_length)
with metric_col_4:
    st.metric("Identity", f"{selected_result.percent_identity:.1f}%")
with metric_col_5:
    st.metric("Conflict columns", selected_result.conflict_count)

if selected_result.aligned_left:
    st.subheader("Gapped assembly")
    st.code(format_assembly_block(selected_result), wrap_lines=False)

conflict_rows = assembly_conflicts_to_rows(selected_result)
selected_conflict_row: dict[str, object] | None = None
if conflict_rows:
    st.subheader("Consensus support / conflict columns")
    st.dataframe(conflict_rows, hide_index=True, width="stretch")
    selected_conflict_index = cast(
        int,
        st.selectbox(
            "Center chromatograms on conflict",
            options=list(range(len(conflict_rows))),
            format_func=lambda index: _conflict_option_label(conflict_rows[index]),
            key=_SELECTED_CONFLICT_WIDGET_KEY,
        ),
    )
    selected_conflict_row = conflict_rows[selected_conflict_index]

left_source_filename, right_source_filename = selected_definition.source_filenames
left_raw_record = prepared_batch.parsed_records[left_source_filename]
right_raw_record = prepared_batch.parsed_records[right_source_filename]
left_trim_result = prepared_batch.trim_results[left_source_filename]
right_trim_result = prepared_batch.trim_results[right_source_filename]

left_view = build_chromatogram_view(left_raw_record, left_trim_result)
right_view = build_chromatogram_view(right_raw_record, right_trim_result)

if left_view.is_renderable and right_view.is_renderable:
    st.subheader("Chromatograms")
    theme_type = str(getattr(getattr(st.context, "theme", None), "type", "light"))
    left_figure = build_chromatogram_figure(left_view, theme_type=theme_type)
    right_figure = build_chromatogram_figure(right_view, theme_type=theme_type)

    left_trace_x = (
        None
        if selected_conflict_row is None
        else selected_conflict_row.get("left_trace_x")
    )
    right_trace_x = (
        None
        if selected_conflict_row is None
        else selected_conflict_row.get("right_trace_x")
    )
    if isinstance(left_trace_x, (int, float)):
        left_figure.update_xaxes(
            range=_centered_x_range(left_view, center=float(left_trace_x))
        )
    if isinstance(right_trace_x, (int, float)):
        right_figure.update_xaxes(
            range=_centered_x_range(right_view, center=float(right_trace_x))
        )

    left_figure.update_layout(height=420, margin={"l": 24, "r": 24, "t": 24, "b": 24})
    right_figure.update_layout(
        height=420,
        margin={"l": 24, "r": 24, "t": 24, "b": 24},
    )

    left_plot_col, right_plot_col = st.columns(2)
    with left_plot_col:
        st.caption(f"Left: {selected_result.left_display_name}")
        st.plotly_chart(
            left_figure,
            width="stretch",
            config={"scrollZoom": False},
        )
    with right_plot_col:
        st.caption(f"Right: {selected_result.right_display_name}")
        st.plotly_chart(
            right_figure,
            width="stretch",
            config={"scrollZoom": False},
        )

consensus_record = selected_computed_assembly.consensus_record
if consensus_record is not None and consensus_record.sequence:
    consensus_fasta = to_fasta(consensus_record, line_width=None)
    st.subheader("Consensus FASTA")
    st.code(consensus_fasta, wrap_lines=True)
    st.download_button(
        "Download consensus FASTA",
        data=consensus_fasta,
        file_name=f"{consensus_record.name}.fasta",
        mime="text/plain",
    )
