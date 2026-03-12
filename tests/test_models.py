from abi_sauce.models import SequenceUpload


def test_sequence_upload_suffix_is_lowercase_without_dot() -> None:
    upload = SequenceUpload(filename="Sample.AB1", content=b"ABC")
    assert upload.suffix == "ab1"


def test_sequence_upload_size_bytes_matches_content_length() -> None:
    upload = SequenceUpload(filename="sample.ab1", content=b"12345")
    assert upload.size_bytes == 5
