# tests/test_trim_state.py
from __future__ import annotations

from abi_sauce.trim_state import (
    DEFAULT_BATCH_TRIM_CONFIG,
    BatchTrimState,
    apply_submitted_trim_config,
    build_record_annotations,
    resolve_active_trim_config,
    resolve_batch_trim_inputs,
)
from abi_sauce.trimming import TrimConfig


def test_apply_submitted_trim_config_removes_selected_noop_override() -> None:
    trim_state = BatchTrimState(
        trim_scope="selected",
        trim_configs_by_record={"a.ab1": TrimConfig(left_trim=2)},
    )

    updated_trim_state = apply_submitted_trim_config(
        trim_state,
        selected_record_name="a.ab1",
        submitted_trim_config=TrimConfig(),
    )

    assert dict(updated_trim_state.trim_configs_by_record) == {}


def test_build_record_annotations_marks_only_overridden_records() -> None:
    annotations = build_record_annotations(
        ["a.ab1", "b.ab1", "c.ab1"],
        {
            "b.ab1": TrimConfig(left_trim=1),
            "c.ab1": TrimConfig(right_trim=1),
        },
    )

    assert annotations.display_labels_by_record == {
        "a.ab1": "a.ab1",
        "b.ab1": "b.ab1 *",
        "c.ab1": "c.ab1 *",
    }
    assert annotations.overridden_record_names == frozenset({"b.ab1", "c.ab1"})
    assert annotations.overridden_count == 2
    assert annotations.custom_trim_flags_by_record == {
        "a.ab1": False,
        "b.ab1": True,
        "c.ab1": True,
    }


def test_build_record_annotations_ignores_stale_override_records() -> None:
    annotations = build_record_annotations(
        ["a.ab1", "b.ab1"],
        {
            "b.ab1": TrimConfig(left_trim=1),
            "stale.ab1": TrimConfig(right_trim=1),
        },
    )

    assert annotations.display_labels_by_record == {
        "a.ab1": "a.ab1",
        "b.ab1": "b.ab1 *",
    }
    assert annotations.overridden_record_names == frozenset({"b.ab1"})
    assert annotations.overridden_count == 1
    assert annotations.custom_trim_flags_by_record == {
        "a.ab1": False,
        "b.ab1": True,
    }


def test_resolve_batch_trim_inputs_uses_only_global_config_in_all_mode() -> None:
    trim_state = BatchTrimState(
        trim_scope="all",
        global_trim_config=TrimConfig(left_trim=3),
        trim_configs_by_record={"a.ab1": TrimConfig(right_trim=2)},
    )

    resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)

    assert resolved_trim_inputs.default_trim_config == TrimConfig(left_trim=3)
    assert resolved_trim_inputs.trim_configs_by_name == {}


def test_resolve_batch_trim_inputs_preserves_batch_default_in_selected_mode() -> None:
    trim_state = BatchTrimState(
        trim_scope="selected",
        global_trim_config=TrimConfig(left_trim=3),
        trim_configs_by_record={"a.ab1": TrimConfig(right_trim=2)},
    )

    resolved_trim_inputs = resolve_batch_trim_inputs(trim_state)

    assert resolved_trim_inputs.default_trim_config == TrimConfig(left_trim=3)
    assert resolved_trim_inputs.trim_configs_by_name == {
        "a.ab1": TrimConfig(right_trim=2),
    }


def test_resolve_active_trim_config_tracks_scope_and_selected_record() -> None:
    all_scope_state = BatchTrimState(
        trim_scope="all",
        global_trim_config=TrimConfig(left_trim=3),
        trim_configs_by_record={"b.ab1": TrimConfig(right_trim=2)},
    )
    selected_scope_state = BatchTrimState(
        trim_scope="selected",
        global_trim_config=all_scope_state.global_trim_config,
        trim_configs_by_record=all_scope_state.trim_configs_by_record,
    )

    assert resolve_active_trim_config(
        selected_scope_state,
        selected_record_name="a.ab1",
    ) == TrimConfig(left_trim=3)
    assert resolve_active_trim_config(
        selected_scope_state,
        selected_record_name="b.ab1",
    ) == TrimConfig(right_trim=2)
    assert resolve_active_trim_config(
        all_scope_state,
        selected_record_name="b.ab1",
    ) == TrimConfig(left_trim=3)


def test_resolve_active_trim_config_uses_builtin_default_when_global_config_is_none() -> (
    None
):
    trim_state = BatchTrimState(
        trim_scope="selected",
        global_trim_config=None,
        trim_configs_by_record={},
    )

    assert (
        resolve_active_trim_config(
            trim_state,
            selected_record_name="a.ab1",
        )
        == DEFAULT_BATCH_TRIM_CONFIG
    )
