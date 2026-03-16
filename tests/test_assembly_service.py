from __future__ import annotations

import io
import json
import zipfile

from abi_sauce.assembly_types import AssemblyConfig, MultiAssemblyResult
from abi_sauce.assembly_state import AssemblyDefinition
from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData
from abi_sauce.services.assembly_compute import (
    compute_saved_assemblies,
    compute_saved_assembly,
)
from abi_sauce.services.assembly_export import (
    accepted_consensus_records,
    prepare_assembly_download,
    select_assembly_export,
)
from abi_sauce.services.batch_parse import ParsedBatch, build_batch_signature
from abi_sauce.services.batch_trim import apply_trim_configs
from abi_sauce.trimming import TrimConfig


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    trace_length = 200
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic assembly record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=TraceData(
            channels={
                "DATA9": [1] * trace_length,
                "DATA10": [1] * trace_length,
                "DATA11": [1] * trace_length,
                "DATA12": [1] * trace_length,
            },
            base_positions=base_positions
            or list(range(5, 5 + (10 * len(sequence)), 10)),
            channel_order="GATC",
        ),
    )


def make_prepared_batch():
    uploads = (
        SequenceUpload(filename="left.ab1", content=b"left"),
        SequenceUpload(filename="right.ab1", content=b"right"),
        SequenceUpload(filename="short.ab1", content=b"short"),
    )
    parsed_batch = ParsedBatch(
        uploads=uploads,
        parsed_records={
            "left.ab1": make_record(
                name="left",
                sequence="CCCCAAAAC",
                qualities=[40] * 9,
            ),
            "right.ab1": make_record(
                name="right",
                sequence="GTTTT",
                qualities=[40] * 5,
            ),
            "short.ab1": make_record(
                name="short",
                sequence="CCCCGGGG",
                qualities=[40] * 8,
            ),
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )
    return apply_trim_configs(
        parsed_batch,
        trim_configs_by_name={"short.ab1": TrimConfig()},
    )


def test_compute_saved_assembly_returns_named_consensus_record() -> None:
    prepared_batch = make_prepared_batch()
    definition = AssemblyDefinition(
        assembly_id="assembly-good",
        name="Amplicon A",
        source_filenames=("left.ab1", "right.ab1"),
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
    )

    computed_assembly = compute_saved_assembly(prepared_batch, definition)

    assert computed_assembly.status == "ok"
    assert computed_assembly.result is not None
    assert computed_assembly.result.accepted is True
    assert computed_assembly.consensus_record is not None
    assert computed_assembly.consensus_record.name == "Amplicon A"
    assert computed_assembly.consensus_record.sequence == "CCCCAAAAC"


def test_compute_saved_assembly_returns_multi_result_for_multi_definition() -> None:
    prepared_batch = make_prepared_batch()
    definition = AssemblyDefinition(
        assembly_id="assembly-multi",
        name="Multi",
        source_filenames=("left.ab1", "right.ab1", "short.ab1"),
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        engine_kind="multi",
    )

    computed_assembly = compute_saved_assembly(prepared_batch, definition)

    assert computed_assembly.status == "ok"
    assert isinstance(computed_assembly.result, MultiAssemblyResult)
    assert computed_assembly.consensus_record is not None
    assert computed_assembly.consensus_record.sequence == "CCCCAAAAC"
    assert computed_assembly.result.included_member_count == 2
    assert computed_assembly.result.excluded_member_count == 1
    assert computed_assembly.result.members[2].included is False


def test_compute_saved_assembly_rejects_invalid_pairwise_definition() -> None:
    prepared_batch = make_prepared_batch()
    definition = AssemblyDefinition(
        assembly_id="assembly-invalid",
        name="Too many reads",
        source_filenames=("left.ab1", "right.ab1", "short.ab1"),
        config=AssemblyConfig(),
    )

    computed_assembly = compute_saved_assembly(prepared_batch, definition)

    assert computed_assembly.status == "invalid_definition"
    assert computed_assembly.result is None
    assert computed_assembly.consensus_record is None
    assert computed_assembly.status_reason is not None
    assert "exactly 2 reads" in computed_assembly.status_reason


def test_compute_saved_assemblies_and_accepted_consensus_records_preserve_order() -> (
    None
):
    prepared_batch = make_prepared_batch()
    definitions = (
        AssemblyDefinition(
            assembly_id="assembly-good",
            name="Amplicon A",
            source_filenames=("left.ab1", "right.ab1"),
            config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        ),
        AssemblyDefinition(
            assembly_id="assembly-rejected",
            name="Rejected",
            source_filenames=("left.ab1", "short.ab1"),
            config=AssemblyConfig(min_overlap_length=9, min_percent_identity=95.0),
        ),
    )

    computed_assemblies = compute_saved_assemblies(prepared_batch, definitions)
    consensus_records = accepted_consensus_records(computed_assemblies)

    assert tuple(computed_assemblies) == ("assembly-good", "assembly-rejected")
    assert tuple(record.name for record in consensus_records) == ("Amplicon A",)


def test_select_assembly_export_can_include_rejected_consensus_records() -> None:
    prepared_batch = make_prepared_batch()
    definitions = (
        AssemblyDefinition(
            assembly_id="assembly-good",
            name="Amplicon A",
            source_filenames=("left.ab1", "right.ab1"),
            config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        ),
        AssemblyDefinition(
            assembly_id="assembly-rejected",
            name="Rejected",
            source_filenames=("left.ab1", "short.ab1"),
            config=AssemblyConfig(min_overlap_length=9, min_percent_identity=95.0),
        ),
    )
    computed_assemblies = compute_saved_assemblies(prepared_batch, definitions)

    accepted_only = select_assembly_export(
        computed_assemblies,
        selected_ids=("assembly-good", "assembly-rejected"),
        require_accepted=True,
    )
    include_rejected = select_assembly_export(
        computed_assemblies,
        selected_ids=("assembly-good", "assembly-rejected"),
        require_accepted=False,
    )

    assert tuple(record.name for record in accepted_only.eligible_records) == (
        "Amplicon A",
    )
    assert tuple(record.name for record in include_rejected.eligible_records) == (
        "Amplicon A",
        "Rejected",
    )
    assert accepted_only.ineligible_reasons == (
        ("Rejected", ("overlap length below threshold (1 < 9)",)),
    )


def test_prepare_assembly_download_builds_concatenated_fasta() -> None:
    prepared_batch = make_prepared_batch()
    computed_assemblies = compute_saved_assemblies(
        prepared_batch,
        (
            AssemblyDefinition(
                assembly_id="assembly-good",
                name="Amplicon A",
                source_filenames=("left.ab1", "right.ab1"),
                config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
            ),
        ),
    )

    artifact = prepare_assembly_download(
        computed_assemblies,
        selected_ids=("assembly-good",),
        concatenate_batch=True,
        filename_stem="consensus-batch",
    )

    assert artifact.is_downloadable is True
    assert artifact.filename == "consensus-batch.fasta"
    assert artifact.mime == "text/plain"
    assert artifact.data == ">Amplicon A\nCCCCAAAAC\n"


def test_prepare_assembly_download_builds_zip_manifest() -> None:
    prepared_batch = make_prepared_batch()
    computed_assemblies = compute_saved_assemblies(
        prepared_batch,
        (
            AssemblyDefinition(
                assembly_id="assembly-good",
                name="Amplicon A",
                source_filenames=("left.ab1", "right.ab1"),
                config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
            ),
            AssemblyDefinition(
                assembly_id="assembly-rejected",
                name="Rejected",
                source_filenames=("left.ab1", "short.ab1"),
                config=AssemblyConfig(min_overlap_length=9, min_percent_identity=95.0),
            ),
        ),
    )

    artifact = prepare_assembly_download(
        computed_assemblies,
        selected_ids=("assembly-good", "assembly-rejected"),
        concatenate_batch=False,
        filename_stem="consensus-batch",
        include_manifest=True,
        require_accepted=True,
    )

    assert artifact.is_downloadable is True
    assert artifact.filename == "consensus-batch.zip"
    assert artifact.mime == "application/zip"

    assert isinstance(artifact.data, bytes)
    with zipfile.ZipFile(io.BytesIO(artifact.data)) as zip_file:
        assert zip_file.namelist() == [
            "001_Amplicon_A.fasta",
            "manifest.json",
        ]
        manifest = json.loads(zip_file.read("manifest.json").decode("utf-8"))

    assert manifest["exported_assemblies"][0]["assembly_name"] == "Amplicon A"
    assert manifest["excluded_assemblies"] == [
        {
            "assembly_name": "Rejected",
            "reasons": ["overlap length below threshold (1 < 9)"],
        }
    ]
