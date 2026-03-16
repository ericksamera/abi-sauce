from __future__ import annotations

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.assembly_state import (
    AssemblyDefinition,
    AssemblySessionState,
    clear_assembly_session_state,
    create_assembly_definition,
    delete_assembly_definition,
    read_assembly_session_state,
    set_export_selected_assembly_ids,
    set_selected_assembly_id,
    suggest_assembly_name,
    sync_assembly_session_state,
    update_assembly_definition,
)
from abi_sauce.models import SequenceUpload
from abi_sauce.services.batch_parse import build_batch_signature


def make_signature(*filenames: str):
    uploads = tuple(
        SequenceUpload(filename=filename, content=filename.encode("utf-8"))
        for filename in filenames
    )
    return build_batch_signature(uploads)


def test_read_assembly_session_state_defaults_to_empty() -> None:
    session_state: dict[str, object] = {}

    assembly_state = read_assembly_session_state(session_state)

    assert assembly_state.batch_signature is None
    assert assembly_state.assemblies_by_id == {}
    assert assembly_state.selected_assembly_id is None
    assert assembly_state.export_selected_assembly_ids == frozenset()


def test_create_update_and_delete_assembly_definition_round_trip() -> None:
    session_state: dict[str, object] = {}
    initial_state = AssemblySessionState(
        batch_signature=make_signature("a.ab1", "b.ab1")
    )
    session_state["abi_sauce.assembly.batch_signature"] = initial_state.batch_signature

    created = create_assembly_definition(
        session_state,
        name="My assembly",
        source_filenames=("a.ab1", "b.ab1"),
        config=AssemblyConfig(min_overlap_length=40),
    )
    created_state = read_assembly_session_state(session_state)

    assert created_state.selected_assembly_id == created.assembly_id
    assert created_state.assemblies_by_id[created.assembly_id].name == "My assembly"
    assert (
        created_state.assemblies_by_id[created.assembly_id].config.min_overlap_length
        == 40
    )

    updated = update_assembly_definition(
        session_state,
        assembly_id=created.assembly_id,
        name="Edited assembly",
        source_filenames=("b.ab1", "a.ab1"),
        config=AssemblyConfig(min_overlap_length=30),
    )
    updated_state = read_assembly_session_state(session_state)

    assert updated.assembly_id == created.assembly_id
    assert updated_state.assemblies_by_id[created.assembly_id].name == "Edited assembly"
    assert updated_state.assemblies_by_id[created.assembly_id].source_filenames == (
        "b.ab1",
        "a.ab1",
    )

    delete_assembly_definition(session_state, assembly_id=created.assembly_id)
    deleted_state = read_assembly_session_state(session_state)

    assert deleted_state.assemblies_by_id == {}
    assert deleted_state.selected_assembly_id is None


def test_sync_assembly_session_state_prunes_stale_definitions_and_export_selection() -> (
    None
):
    batch_signature = make_signature("a.ab1", "b.ab1")
    stale_definition = AssemblyDefinition(
        assembly_id="assembly-stale",
        name="Stale",
        source_filenames=("stale.ab1", "b.ab1"),
        config=AssemblyConfig(),
    )
    kept_definition = AssemblyDefinition(
        assembly_id="assembly-keep",
        name="Keep",
        source_filenames=("a.ab1", "b.ab1"),
        config=AssemblyConfig(),
    )
    session_state: dict[str, object] = {
        "abi_sauce.assembly.batch_signature": batch_signature,
        "abi_sauce.assembly.assemblies_by_id": {
            stale_definition.assembly_id: stale_definition,
            kept_definition.assembly_id: kept_definition,
        },
        "abi_sauce.assembly.selected_id": stale_definition.assembly_id,
        "abi_sauce.assembly.export_selected_ids": (
            stale_definition.assembly_id,
            kept_definition.assembly_id,
        ),
    }

    assembly_state = sync_assembly_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("a.ab1", "b.ab1"),
    )

    assert tuple(assembly_state.assemblies_by_id) == (kept_definition.assembly_id,)
    assert assembly_state.selected_assembly_id == kept_definition.assembly_id
    assert assembly_state.export_selected_assembly_ids == frozenset(
        {kept_definition.assembly_id}
    )


def test_setters_and_clear_assembly_session_state() -> None:
    session_state: dict[str, object] = {}
    created = create_assembly_definition(
        session_state,
        name="",
        source_filenames=("a.ab1", "b.ab1"),
        config=AssemblyConfig(),
    )

    set_selected_assembly_id(session_state, created.assembly_id)
    set_export_selected_assembly_ids(session_state, [created.assembly_id])

    assembly_state = read_assembly_session_state(session_state)
    assert assembly_state.selected_assembly_id == created.assembly_id
    assert assembly_state.export_selected_assembly_ids == frozenset(
        {created.assembly_id}
    )

    clear_assembly_session_state(session_state)
    assert read_assembly_session_state(session_state) == AssemblySessionState()


def test_create_assembly_definition_preserves_multi_engine_kind() -> None:
    session_state: dict[str, object] = {}

    created = create_assembly_definition(
        session_state,
        name="Multi assembly",
        source_filenames=("a.ab1", "b.ab1", "c.ab1"),
        config=AssemblyConfig(),
        engine_kind="multi",
    )

    assembly_state = read_assembly_session_state(session_state)

    assert created.engine_kind == "multi"
    assert assembly_state.assemblies_by_id[created.assembly_id].engine_kind == "multi"


def test_suggest_assembly_name_dedupes_against_existing_names() -> None:
    assert (
        suggest_assembly_name(
            ("sample_f.ab1", "sample_r.ab1"),
            existing_names=("sample_f + sample_r",),
        )
        == "sample_f + sample_r (2)"
    )
