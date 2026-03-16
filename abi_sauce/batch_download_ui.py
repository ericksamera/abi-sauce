from __future__ import annotations

from typing import cast

import streamlit as st

from abi_sauce.batch import ExportFormat
from abi_sauce.exceptions import ExportError
from abi_sauce.services.batch_export import prepare_batch_download
from abi_sauce.services.batch_trim import PreparedBatch


def render_batch_download_controls(
    *,
    prepared_batch: PreparedBatch,
    key_prefix: str,
    default_filename_stem: str = "",
    button_label: str = "Download trimmed batch",
    compact: bool = False,
) -> None:
    """Render reusable batch-download controls for the active trimmed batch."""
    export_format = cast(
        ExportFormat,
        st.selectbox(
            "Format",
            options=["fasta", "fastq"],
            key=f"{key_prefix}.export_format",
        ),
    )

    concatenate_batch = st.checkbox(
        "Concatenate entries into a single file",
        value=True,
        key=f"{key_prefix}.concatenate_batch",
    )
    if not concatenate_batch:
        st.caption("Individual records will be bundled into a ZIP archive.")

    filename_stem = st.text_input(
        "Output filename stem",
        value=default_filename_stem,
        key=f"{key_prefix}.filename_stem",
    )
    exclude_failed_min_length = st.checkbox(
        "Exclude sequences failing minimum length",
        key=f"{key_prefix}.exclude_failed_min_length",
        value=True,
    )

    try:
        download_artifact = prepare_batch_download(
            prepared_batch,
            export_format=export_format,
            concatenate_batch=bool(concatenate_batch),
            filename_stem=filename_stem,
            require_min_length=bool(exclude_failed_min_length),
        )
    except ExportError as exc:
        st.warning(str(exc))
        return

    if compact:
        st.caption(
            f"{len(download_artifact.eligible_records)} exportable record(s) · "
            f"{len(download_artifact.ineligible_reasons)} excluded"
        )
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
        return

    @st.fragment
    def _download_fragment() -> None:
        st.download_button(
            label=button_label,
            data=download_artifact.data,
            file_name=download_artifact.filename,
            mime=download_artifact.mime,
            key=f"{key_prefix}.download_button",
        )

    _download_fragment()
