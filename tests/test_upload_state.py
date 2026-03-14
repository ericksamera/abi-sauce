from __future__ import annotations

import pytest

from abi_sauce.models import SequenceRecord, SequenceUpload
from abi_sauce.services.batch import ParsedBatch, build_batch_signature
from abi_sauce.upload_state import (
    clear_active_batch,
    get_active_batch_signature,
    get_active_parsed_batch,
    get_active_uploads,
    merge_uploads,
    read_active_batch_state,
    set_active_parsed_batch,
    set_active_uploads,
    update_active_parsed_record,
)


def make_upload(filename: str, content: bytes = b"fake") -> SequenceUpload:
    return SequenceUpload(filename=filename, content=content)


def make_record(name: str) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="test record",
        sequence="ACGT",
        source_format="abi",
        qualities=[40, 41, 42, 43],
    )


def make_parsed_batch() -> ParsedBatch:
    uploads = (
        make_upload("a.ab1", b"aa"),
        make_upload("b.ab1", b"bbb"),
    )
    return ParsedBatch(
        uploads=uploads,
        parsed_records={
            "a.ab1": make_record("trace_a"),
            "b.ab1": make_record("trace_b"),
        },
        parse_errors={},
        signature=build_batch_signature(uploads),
    )


def test_active_batch_state_defaults_to_empty() -> None:
    session_state: dict[str, object] = {}

    active_batch_state = read_active_batch_state(session_state)

    assert active_batch_state.uploads == ()
    assert active_batch_state.parsed_batch is None
    assert active_batch_state.signature is None
    assert get_active_uploads(session_state) == ()
    assert get_active_parsed_batch(session_state) is None
    assert get_active_batch_signature(session_state) is None


def test_set_active_uploads_stores_normalized_uploads_without_parsed_batch() -> None:
    session_state: dict[str, object] = {}
    uploads = (
        make_upload("a.ab1", b"aa"),
        make_upload("b.ab1", b"bbb"),
    )

    set_active_uploads(session_state, uploads)

    assert get_active_uploads(session_state) == uploads
    assert get_active_parsed_batch(session_state) is None
    assert get_active_batch_signature(session_state) == build_batch_signature(uploads)


def test_set_active_parsed_batch_persists_uploads_and_batch_together() -> None:
    session_state: dict[str, object] = {}
    parsed_batch = make_parsed_batch()

    set_active_parsed_batch(session_state, parsed_batch)

    active_batch_state = read_active_batch_state(session_state)

    assert active_batch_state.uploads == parsed_batch.uploads
    assert active_batch_state.parsed_batch == parsed_batch
    assert active_batch_state.signature == parsed_batch.signature
    assert get_active_uploads(session_state) == parsed_batch.uploads
    assert get_active_parsed_batch(session_state) == parsed_batch
    assert get_active_batch_signature(session_state) == parsed_batch.signature


def test_get_active_uploads_prefers_parsed_batch_when_stored_uploads_are_stale() -> (
    None
):
    parsed_batch = make_parsed_batch()
    session_state: dict[str, object] = {
        "abi_sauce.active_uploads": (make_upload("stale.ab1", b"xxxx"),),
        "abi_sauce.active_parsed_batch": parsed_batch,
    }

    assert get_active_uploads(session_state) == parsed_batch.uploads
    assert get_active_batch_signature(session_state) == parsed_batch.signature


def test_set_active_uploads_rejects_mismatched_parsed_batch() -> None:
    session_state: dict[str, object] = {}
    parsed_batch = make_parsed_batch()

    with pytest.raises(ValueError, match="must match uploads"):
        set_active_uploads(
            session_state,
            (make_upload("different.ab1", b"zz"),),
            parsed_batch=parsed_batch,
        )


def test_update_active_parsed_record_replaces_one_record_in_session_batch() -> None:
    session_state: dict[str, object] = {}
    parsed_batch = make_parsed_batch()
    set_active_parsed_batch(session_state, parsed_batch)

    replacement_record = make_record("trace_a_rc")
    replacement_record.orientation = "reverse_complement"

    updated_batch = update_active_parsed_record(
        session_state,
        source_filename="a.ab1",
        record=replacement_record,
    )

    assert updated_batch.uploads == parsed_batch.uploads
    assert updated_batch.parse_errors == parsed_batch.parse_errors
    assert updated_batch.signature == parsed_batch.signature
    assert get_active_parsed_batch(session_state) == updated_batch
    assert updated_batch.parsed_records["a.ab1"].orientation == "reverse_complement"
    assert updated_batch.parsed_records["b.ab1"] == parsed_batch.parsed_records["b.ab1"]


def test_update_active_parsed_record_requires_an_active_parsed_batch() -> None:
    session_state: dict[str, object] = {}

    with pytest.raises(ValueError, match="No active parsed batch is available."):
        update_active_parsed_record(
            session_state,
            source_filename="a.ab1",
            record=make_record("trace_a"),
        )


def test_merge_uploads_appends_new_files_and_replaces_name_collisions() -> None:
    merged_uploads = merge_uploads(
        (
            make_upload("a.ab1", b"old-a"),
            make_upload("b.ab1", b"old-b"),
        ),
        (
            make_upload("b.ab1", b"new-b"),
            make_upload("c.ab1", b"new-c"),
        ),
    )

    assert tuple(upload.filename for upload in merged_uploads) == (
        "a.ab1",
        "b.ab1",
        "c.ab1",
    )
    assert merged_uploads[0].content == b"old-a"
    assert merged_uploads[1].content == b"new-b"
    assert merged_uploads[2].content == b"new-c"


def test_clear_active_batch_removes_all_session_entries() -> None:
    session_state: dict[str, object] = {}
    parsed_batch = make_parsed_batch()
    set_active_parsed_batch(session_state, parsed_batch)

    clear_active_batch(session_state)

    assert read_active_batch_state(session_state).uploads == ()
    assert get_active_parsed_batch(session_state) is None
    assert get_active_batch_signature(session_state) is None
