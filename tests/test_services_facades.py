from __future__ import annotations

from abi_sauce.exceptions import AbiParseError
from abi_sauce.parsers.abi import parse_ab1_upload
from abi_sauce.services import alignment as alignment_service
from abi_sauce.services import assembly as assembly_service
from abi_sauce.services import batch as batch_service
from abi_sauce.services import (
    alignment_compute,
    assembly_compute,
    assembly_export,
    batch_export,
    batch_parse,
    batch_trim,
)


def test_assembly_service_facade_reexports_dedicated_modules() -> None:
    assert (
        assembly_service.compute_saved_assemblies
        is assembly_compute.compute_saved_assemblies
    )
    assert (
        assembly_service.compute_saved_assembly
        is assembly_compute.compute_saved_assembly
    )
    assert (
        assembly_service.compute_saved_multi_assembly
        is assembly_compute.compute_saved_multi_assembly
    )
    assert (
        assembly_service.accepted_consensus_records
        is assembly_export.accepted_consensus_records
    )
    assert (
        assembly_service.prepare_assembly_download
        is assembly_export.prepare_assembly_download
    )
    assert (
        assembly_service.select_assembly_export
        is assembly_export.select_assembly_export
    )


def test_batch_service_facade_reexports_dedicated_modules() -> None:
    assert batch_service.AbiParseError is AbiParseError
    assert batch_service.parse_ab1_upload is parse_ab1_upload
    assert batch_service.ParsedBatch is batch_parse.ParsedBatch
    assert batch_service.PreparedBatch is batch_trim.PreparedBatch
    assert batch_service.BatchDownloadArtifact is batch_export.BatchDownloadArtifact
    assert batch_service.BatchExportSelection is batch_export.BatchExportSelection
    assert batch_service.build_batch_signature is batch_parse.build_batch_signature
    assert (
        batch_service.normalize_uploaded_files is batch_parse.normalize_uploaded_files
    )
    assert (
        batch_service.replace_parsed_batch_record
        is batch_parse.replace_parsed_batch_record
    )
    assert batch_service.apply_trim_config is batch_trim.apply_trim_config
    assert batch_service.apply_trim_configs is batch_trim.apply_trim_configs
    assert batch_service.prepare_batch_download is batch_export.prepare_batch_download
    assert batch_service.select_batch_export is batch_export.select_batch_export


def test_alignment_service_facade_reexports_dedicated_modules() -> None:
    assert (
        alignment_service.compute_saved_alignment
        is alignment_compute.compute_saved_alignment
    )
    assert (
        alignment_service.compute_saved_alignments
        is alignment_compute.compute_saved_alignments
    )
