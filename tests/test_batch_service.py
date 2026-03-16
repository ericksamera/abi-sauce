from __future__ import annotations

import io
import zipfile

import hashlib
import pytest

from abi_sauce.exceptions import AbiParseError, ExportError
from abi_sauce.models import SequenceOrientation, SequenceRecord, SequenceUpload
from abi_sauce.services.batch_export import prepare_batch_download, select_batch_export
from abi_sauce.services.batch_parse import (
    ParsedBatch,
    build_batch_signature,
    normalize_uploaded_files,
    parse_uploaded_batch,
    replace_parsed_batch_record,
)
from abi_sauce.services.batch_trim import apply_trim_config, apply_trim_configs
from abi_sauce.trimming import TrimConfig


def expected_digest(content: bytes) -> str:
    return hashlib.blake2b(content, digest_size=16).hexdigest()


class FakeUploadedFile:
    def __init__(self, name: str, content: bytes = b"fake") -> None:
        self.name = name
        self._content = content

    def getvalue(self) -> bytes:
        return self._content


def make_record(
    name: str,
    *,
    sequence: str = "ACGT",
    qualities: list[int] | None = None,
    orientation: SequenceOrientation = "forward",
) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="test record",
        sequence=sequence,
        source_format="abi",
        orientation=orientation,
        qualities=qualities,
    )


def test_normalize_uploaded_files_sorts_and_wraps() -> None:
    uploads = normalize_uploaded_files(
        [
            FakeUploadedFile("b.ab1", b"bbb"),
            FakeUploadedFile("a.ab1", b"aa"),
        ]
    )

    assert tuple(upload.filename for upload in uploads) == ("a.ab1", "b.ab1")
    assert tuple(upload.size_bytes for upload in uploads) == (2, 3)


def test_parse_uploaded_batch_collects_records_and_errors() -> None:
    def fake_parse_ab1_upload(upload):
        if upload.filename == "broken.ab1":
            raise AbiParseError("Failed to parse ABI file: broken.ab1")
        return make_record(upload.filename.removesuffix(".ab1"))

    parsed_batch = parse_uploaded_batch(
        [
            FakeUploadedFile("b.ab1", b"bbb"),
            FakeUploadedFile("broken.ab1", b"x"),
            FakeUploadedFile("a.ab1", b"aa"),
        ],
        parse_upload=fake_parse_ab1_upload,
    )

    assert tuple(upload.filename for upload in parsed_batch.uploads) == (
        "a.ab1",
        "b.ab1",
        "broken.ab1",
    )
    assert tuple(parsed_batch.parsed_records) == ("a.ab1", "b.ab1")
    assert parsed_batch.parse_errors == {
        "broken.ab1": "Failed to parse ABI file: broken.ab1",
    }
    assert parsed_batch.signature == (
        ("a.ab1", 2, expected_digest(b"aa")),
        ("b.ab1", 3, expected_digest(b"bbb")),
        ("broken.ab1", 1, expected_digest(b"x")),
    )


def test_build_batch_signature_changes_when_content_changes_but_name_and_size_do_not() -> (
    None
):
    first = (SequenceUpload(filename="trace.ab1", content=b"AAAA"),)
    second = (SequenceUpload(filename="trace.ab1", content=b"TTTT"),)

    assert first[0].filename == second[0].filename
    assert first[0].size_bytes == second[0].size_bytes
    assert build_batch_signature(first) != build_batch_signature(second)


def make_parsed_batch() -> ParsedBatch:
    uploads = (
        SequenceUpload(filename="a.ab1", content=b"aaaaaa"),
        SequenceUpload(filename="b.ab1", content=b"bbbb"),
        SequenceUpload(filename="broken.ab1", content=b"x"),
    )
    return ParsedBatch(
        uploads=uploads,
        parsed_records={
            "a.ab1": make_record(
                "trace_a",
                sequence="ACGTAC",
                qualities=[40, 41, 42, 43, 44, 45],
            ),
            "b.ab1": make_record(
                "trace_b",
                sequence="ACGT",
                qualities=None,
            ),
        },
        parse_errors={"broken.ab1": "Failed to parse ABI file: broken.ab1"},
        signature=build_batch_signature(uploads),
    )


def test_replace_parsed_batch_record_replaces_one_record_and_preserves_batch_metadata() -> (
    None
):
    parsed_batch = make_parsed_batch()
    replacement_record = make_record(
        "trace_a_rc",
        sequence="TTTT",
        qualities=[30, 31, 32, 33],
        orientation="reverse_complement",
    )

    updated_batch = replace_parsed_batch_record(
        parsed_batch,
        source_filename="a.ab1",
        record=replacement_record,
    )

    assert updated_batch.uploads == parsed_batch.uploads
    assert updated_batch.parse_errors == parsed_batch.parse_errors
    assert updated_batch.signature == parsed_batch.signature
    assert updated_batch.parsed_records["a.ab1"] == replacement_record
    assert updated_batch.parsed_records["b.ab1"] == parsed_batch.parsed_records["b.ab1"]
    assert parsed_batch.parsed_records["a.ab1"].name == "trace_a"


def test_apply_trim_configs_supports_per_record_configs() -> None:
    prepared_batch = apply_trim_configs(
        make_parsed_batch(),
        trim_configs_by_name={
            "a.ab1": TrimConfig(left_trim=2),
        },
    )

    assert prepared_batch.trim_results["a.ab1"].record.sequence == "GTAC"
    assert prepared_batch.trim_results["b.ab1"].record.sequence == "ACGT"


def test_apply_trim_configs_supports_default_config_with_per_record_override() -> None:
    prepared_batch = apply_trim_configs(
        make_parsed_batch(),
        default_trim_config=TrimConfig(left_trim=1),
        trim_configs_by_name={
            "b.ab1": TrimConfig(right_trim=1, min_length=4),
        },
    )

    assert prepared_batch.trim_results["a.ab1"].record.sequence == "CGTAC"
    assert prepared_batch.trim_results["b.ab1"].record.sequence == "ACG"
    assert prepared_batch.trim_results["a.ab1"].passed_min_length is True
    assert prepared_batch.trim_results["b.ab1"].passed_min_length is False


def test_apply_trim_config_builds_prepared_batch() -> None:
    parsed_batch = make_parsed_batch()

    prepared_batch = apply_trim_config(
        parsed_batch,
        TrimConfig(left_trim=1, right_trim=1, min_length=3),
    )

    assert prepared_batch.uploads == parsed_batch.uploads
    assert prepared_batch.parsed_records == parsed_batch.parsed_records
    assert prepared_batch.parse_errors == parsed_batch.parse_errors
    assert prepared_batch.signature == parsed_batch.signature
    assert tuple(prepared_batch.trim_results) == ("a.ab1", "b.ab1")
    assert prepared_batch.batch_summary.total_uploaded_files == 3
    assert prepared_batch.batch_summary.failed_files == 1
    assert prepared_batch.batch_summary.trimmed_records == 2
    assert prepared_batch.batch_summary.records_passing_min_length == 1
    assert prepared_batch.batch_summary.records_failing_min_length == 1
    assert (
        prepared_batch.batch_export_policy.fastq_records[0].source_filename == "a.ab1"
    )


def test_select_batch_export_returns_filtered_records_and_reasons() -> None:
    prepared_batch = apply_trim_config(
        make_parsed_batch(),
        TrimConfig(left_trim=1, right_trim=1, min_length=3),
    )

    export_selection = select_batch_export(
        prepared_batch,
        export_format="fastq",
        require_min_length=True,
    )

    assert tuple(record.name for record in export_selection.eligible_records) == (
        "trace_a",
    )
    assert export_selection.ineligible_reasons == (
        ("b.ab1", ("missing per-base qualities", "below minimum length")),
    )


def test_prepare_batch_download_builds_fastq_artifact() -> None:
    prepared_batch = apply_trim_config(
        make_parsed_batch(),
        TrimConfig(left_trim=1, right_trim=1, min_length=3),
    )

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fastq",
        concatenate_batch=True,
        filename_stem="   ",
        require_min_length=True,
    )

    assert artifact.is_downloadable is True
    assert tuple(record.name for record in artifact.eligible_records) == ("trace_a",)
    assert artifact.ineligible_reasons == (
        ("b.ab1", ("missing per-base qualities", "below minimum length")),
    )
    assert artifact.filename == "abi-sauce-batch.fastq"
    assert artifact.mime == "text/plain"
    assert artifact.data == "@trace_a\nCGTA\n+\nJKLM\n"


def test_prepare_batch_download_supports_unwrapped_fasta_output() -> None:
    prepared_batch = apply_trim_config(
        make_parsed_batch(),
        TrimConfig(min_length=1),
    )

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fasta",
        concatenate_batch=True,
        filename_stem="blast-batch",
        require_min_length=True,
        fasta_line_width=None,
    )

    assert artifact.is_downloadable is True
    assert artifact.filename == "blast-batch.fasta"
    assert artifact.mime == "text/plain"
    assert artifact.data == ">trace_a\nACGTAC\n>trace_b\nACGT\n"


def test_prepare_batch_download_applies_orientation_to_fasta_exports() -> None:
    uploads = (SequenceUpload(filename="rc.ab1", content=b"rc"),)
    parsed_batch = ParsedBatch(
        uploads=uploads,
        parsed_records={
            "rc.ab1": make_record(
                "trace_rc",
                sequence="AAGTC",
                orientation="reverse_complement",
            )
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )
    prepared_batch = apply_trim_config(parsed_batch, TrimConfig())

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fasta",
        concatenate_batch=True,
        filename_stem="blast-batch",
        require_min_length=True,
        fasta_line_width=None,
    )

    assert artifact.is_downloadable is True
    assert artifact.data == ">trace_rc\nGACTT\n"


def test_prepare_batch_download_builds_zip_artifact() -> None:
    prepared_batch = apply_trim_config(
        make_parsed_batch(),
        TrimConfig(left_trim=1, right_trim=1, min_length=3),
    )

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fasta",
        concatenate_batch=False,
        filename_stem="trimmed-batch",
        require_min_length=False,
    )

    assert artifact.is_downloadable is True
    assert tuple(record.name for record in artifact.eligible_records) == (
        "trace_a",
        "trace_b",
    )
    assert artifact.ineligible_reasons == ()
    assert artifact.filename == "trimmed-batch.zip"
    assert artifact.mime == "application/zip"
    assert isinstance(artifact.data, bytes)

    with zipfile.ZipFile(io.BytesIO(artifact.data)) as zip_file:
        assert zip_file.namelist() == [
            "001_trace_a.fasta",
            "002_trace_b.fasta",
        ]
        assert zip_file.read("001_trace_a.fasta").decode() == ">trace_a\nCGTA\n"
        assert zip_file.read("002_trace_b.fasta").decode() == ">trace_b\nCG\n"


def test_prepare_batch_download_returns_empty_artifact_when_nothing_is_eligible() -> (
    None
):
    prepared_batch = apply_trim_config(
        make_parsed_batch(),
        TrimConfig(left_trim=1, right_trim=1, min_length=10),
    )

    artifact = prepare_batch_download(
        prepared_batch,
        export_format="fastq",
        concatenate_batch=True,
        filename_stem="trimmed-batch",
        require_min_length=True,
    )

    assert artifact.is_downloadable is False
    assert artifact.eligible_records == ()
    assert artifact.ineligible_reasons == (
        ("a.ab1", ("below minimum length",)),
        ("b.ab1", ("missing per-base qualities", "below minimum length")),
    )
    assert artifact.data == ""
    assert artifact.filename == ""
    assert artifact.mime == ""


def test_prepare_batch_download_propagates_serializer_errors() -> None:
    uploads = (SequenceUpload(filename="bad.ab1", content=b"bad"),)
    parsed_batch = ParsedBatch(
        uploads=uploads,
        parsed_records={
            "bad.ab1": make_record(
                "trace_bad",
                sequence="ACGT",
                qualities=[40, 41, 42, 120],
            )
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )
    prepared_batch = apply_trim_config(parsed_batch, TrimConfig())

    with pytest.raises(
        ExportError,
        match=r"FASTQ export requires PHRED scores between 0 and 93",
    ):
        prepare_batch_download(
            prepared_batch,
            export_format="fastq",
            concatenate_batch=True,
            filename_stem="bad-batch",
        )
