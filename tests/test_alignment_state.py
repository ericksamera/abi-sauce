from __future__ import annotations

from abi_sauce.alignment_state import (
    AlignmentDefinition,
    AlignmentSessionState,
    clear_alignment_session_state,
    create_alignment_definition,
    delete_alignment_definition,
    read_alignment_session_state,
    set_selected_alignment_id,
    suggest_alignment_name,
    sync_alignment_session_state,
    update_alignment_definition,
)
from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceUpload
from abi_sauce.services.batch_parse import build_batch_signature


def make_signature(*filenames: str):
    uploads = tuple(
        SequenceUpload(filename=filename, content=filename.encode("utf-8"))
        for filename in filenames
    )
    return build_batch_signature(uploads)


def test_read_alignment_session_state_defaults_to_empty() -> None:
    session_state: dict[str, object] = {}

    alignment_state = read_alignment_session_state(session_state)

    assert alignment_state.batch_signature is None
    assert alignment_state.alignments_by_id == {}
    assert alignment_state.selected_alignment_id is None


def test_create_update_and_delete_reference_alignment_round_trip() -> None:
    session_state: dict[str, object] = {}
    initial_state = AlignmentSessionState(
        batch_signature=make_signature("trace.ab1"),
    )
    session_state["abi_sauce.alignments.batch_signature"] = (
        initial_state.batch_signature
    )

    created = create_alignment_definition(
        session_state,
        name="Trace vs ref",
        source_filenames=("trace.ab1",),
        engine_kind="reference_single",
        reference_name="amplicon_001",
        reference_text=">amplicon_001\nAACCGGTT\n",
        strand_policy="forward",
    )
    created_state = read_alignment_session_state(session_state)

    assert created_state.selected_alignment_id == created.alignment_id
    assert created_state.alignments_by_id[created.alignment_id].name == "Trace vs ref"
    assert (
        created_state.alignments_by_id[created.alignment_id].reference_name
        == "amplicon_001"
    )

    updated = update_alignment_definition(
        session_state,
        alignment_id=created.alignment_id,
        name="Trace vs ref (edited)",
        source_filenames=("trace.ab1",),
        engine_kind="reference_single",
        reference_name="amplicon_002",
        reference_text=">amplicon_002\nAACCGGTA\n",
        strand_policy="reverse_complement",
    )
    updated_state = read_alignment_session_state(session_state)

    assert updated.alignment_id == created.alignment_id
    assert (
        updated_state.alignments_by_id[created.alignment_id].name
        == "Trace vs ref (edited)"
    )
    assert (
        updated_state.alignments_by_id[created.alignment_id].strand_policy
        == "reverse_complement"
    )

    delete_alignment_definition(session_state, alignment_id=created.alignment_id)
    deleted_state = read_alignment_session_state(session_state)

    assert deleted_state.alignments_by_id == {}
    assert deleted_state.selected_alignment_id is None


def test_sync_alignment_session_state_prunes_stale_definitions() -> None:
    batch_signature = make_signature("trace.ab1", "left.ab1", "right.ab1")
    stale_definition = AlignmentDefinition(
        alignment_id="alignment-stale",
        name="Stale",
        source_filenames=("missing.ab1",),
        engine_kind="reference_single",
        reference_text=">ref\nAACCGGTT\n",
    )
    kept_definition = AlignmentDefinition(
        alignment_id="alignment-keep",
        name="Keep",
        source_filenames=("left.ab1", "right.ab1"),
        engine_kind="pairwise",
        assembly_config=AssemblyConfig(),
    )
    session_state: dict[str, object] = {
        "abi_sauce.alignments.batch_signature": batch_signature,
        "abi_sauce.alignments.alignments_by_id": {
            stale_definition.alignment_id: stale_definition,
            kept_definition.alignment_id: kept_definition,
        },
        "abi_sauce.alignments.selected_id": stale_definition.alignment_id,
    }

    alignment_state = sync_alignment_session_state(
        session_state,
        batch_signature=batch_signature,
        parsed_record_names=("trace.ab1", "left.ab1", "right.ab1"),
    )

    assert tuple(alignment_state.alignments_by_id) == (kept_definition.alignment_id,)
    assert alignment_state.selected_alignment_id == kept_definition.alignment_id


def test_setter_and_clear_alignment_session_state() -> None:
    session_state: dict[str, object] = {}
    created = create_alignment_definition(
        session_state,
        name="",
        source_filenames=("a.ab1", "b.ab1"),
        engine_kind="pairwise",
        assembly_config=AssemblyConfig(),
    )

    set_selected_alignment_id(session_state, created.alignment_id)

    alignment_state = read_alignment_session_state(session_state)
    assert alignment_state.selected_alignment_id == created.alignment_id

    clear_alignment_session_state(session_state)
    cleared_state = read_alignment_session_state(session_state)
    assert cleared_state == AlignmentSessionState()


def test_suggest_alignment_name_derives_pairwise_and_reference_names() -> None:
    pairwise_name = suggest_alignment_name(
        ("left.ab1", "right.ab1"),
        engine_kind="pairwise",
    )
    reference_name = suggest_alignment_name(
        ("trace.ab1",),
        engine_kind="reference_single",
        reference_name="amplicon_001",
    )

    assert pairwise_name == "left + right"
    assert reference_name == "trace vs amplicon_001"
