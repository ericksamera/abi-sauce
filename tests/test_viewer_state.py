from __future__ import annotations

from abi_sauce.services.batch import build_batch_signature
from abi_sauce.trim_state import BatchTrimState
from abi_sauce.trimming import TrimConfig
from abi_sauce.viewer_state import (
    clear_viewer_session_state,
    get_batch_trim_state,
    get_selected_record_name,
    read_viewer_session_state,
    set_batch_trim_state,
    set_selected_record_name,
    sync_viewer_session_state,
)
from abi_sauce.models import SequenceUpload


def make_signature(*filenames: str):
    uploads = tuple(
        SequenceUpload(filename=filename, content=filename.encode("utf-8"))
        for filename in filenames
    )
    return build_batch_signature(uploads)


def test_read_viewer_session_state_defaults_to_empty() -> None:
    session_state: dict[str, object] = {}

    viewer_state = read_viewer_session_state(session_state)

    assert viewer_state.batch_signature is None
    assert viewer_state.selected_record_name is None
    assert viewer_state.trim_state == BatchTrimState()


def test_sync_viewer_session_state_initializes_defaults_for_new_batch() -> None:
    session_state: dict[str, object] = {}
    batch_signature = make_signature("a.ab1", "b.ab1")

    viewer_state = sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )

    assert viewer_state.batch_signature == batch_signature
    assert viewer_state.selected_record_name == "a.ab1"
    assert viewer_state.trim_state == BatchTrimState()
    assert get_selected_record_name(session_state) == "a.ab1"
    assert get_batch_trim_state(session_state) == BatchTrimState()


def test_sync_viewer_session_state_preserves_existing_state_for_same_batch() -> None:
    batch_signature = make_signature("a.ab1", "b.ab1")
    session_state: dict[str, object] = {}

    sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )
    set_selected_record_name(session_state, "b.ab1")
    set_batch_trim_state(
        session_state,
        BatchTrimState(
            trim_scope="selected",
            trim_configs_by_record={"b.ab1": TrimConfig(left_trim=3)},
        ),
    )

    viewer_state = sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )

    assert viewer_state.selected_record_name == "b.ab1"
    assert viewer_state.trim_state.trim_scope == "selected"
    assert viewer_state.trim_state.trim_configs_by_record == {
        "b.ab1": TrimConfig(left_trim=3),
    }


def test_sync_viewer_session_state_prunes_stale_selection_and_overrides() -> None:
    batch_signature = make_signature("a.ab1", "b.ab1")
    session_state: dict[str, object] = {}
    sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )
    set_selected_record_name(session_state, "stale.ab1")
    set_batch_trim_state(
        session_state,
        BatchTrimState(
            trim_scope="selected",
            trim_configs_by_record={
                "b.ab1": TrimConfig(right_trim=2),
                "stale.ab1": TrimConfig(left_trim=1),
            },
        ),
    )

    viewer_state = sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )

    assert viewer_state.selected_record_name == "a.ab1"
    assert viewer_state.trim_state.trim_configs_by_record == {
        "b.ab1": TrimConfig(right_trim=2),
    }


def test_sync_viewer_session_state_resets_trim_state_when_batch_changes() -> None:
    first_signature = make_signature("a.ab1", "b.ab1")
    second_signature = make_signature("c.ab1")
    session_state: dict[str, object] = {}

    sync_viewer_session_state(
        session_state,
        batch_signature=first_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )
    set_selected_record_name(session_state, "b.ab1")
    set_batch_trim_state(
        session_state,
        BatchTrimState(
            trim_scope="selected",
            trim_configs_by_record={"b.ab1": TrimConfig(left_trim=2)},
        ),
    )

    viewer_state = sync_viewer_session_state(
        session_state,
        batch_signature=second_signature,
        parsed_record_names=("c.ab1",),
    )

    assert viewer_state.batch_signature == second_signature
    assert viewer_state.selected_record_name == "c.ab1"
    assert viewer_state.trim_state == BatchTrimState()


def test_clear_viewer_session_state_removes_shared_entries() -> None:
    batch_signature = make_signature("a.ab1")
    session_state: dict[str, object] = {}
    sync_viewer_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1",),
    )
    set_batch_trim_state(
        session_state,
        BatchTrimState(global_trim_config=TrimConfig(left_trim=2)),
    )

    clear_viewer_session_state(session_state)

    assert read_viewer_session_state(session_state).batch_signature is None
    assert get_selected_record_name(session_state) is None
    assert get_batch_trim_state(session_state) == BatchTrimState()
