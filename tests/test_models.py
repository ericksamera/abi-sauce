from abi_sauce.models import SequenceRecord, SequenceUpload


def test_sequence_upload_suffix_is_lowercase_without_dot() -> None:
    upload = SequenceUpload(filename="Sample.AB1", content=b"ABC")
    assert upload.suffix == "ab1"


def test_sequence_upload_size_bytes_matches_content_length() -> None:
    upload = SequenceUpload(filename="sample.ab1", content=b"12345")
    assert upload.size_bytes == 5


def test_sequence_record_orientation_defaults_to_forward() -> None:
    record = SequenceRecord(
        record_id="trace_001",
        name="trace_001",
        description="test",
        sequence="ACGT",
        source_format="abi",
    )

    assert record.orientation == "forward"
