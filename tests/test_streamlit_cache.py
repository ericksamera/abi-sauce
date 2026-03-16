from __future__ import annotations

# ruff: noqa: E402

from dataclasses import replace

import pytest
import streamlit as st

pytest.importorskip("streamlit")

from abi_sauce.assembly_types import AssemblyConfig, MultiAssemblyResult
from abi_sauce.assembly_state import AssemblyDefinition
from abi_sauce.assembly_trace import (
    build_assembly_trace_row_source,
    build_multi_assembly_trace_view,
    build_pairwise_assembly_trace_view,
)
from abi_sauce.models import SequenceRecord, SequenceUpload, TraceData
from abi_sauce.services.assembly_compute import compute_saved_assemblies
from abi_sauce.services.batch_parse import ParsedBatch, build_batch_signature
from abi_sauce.services.batch_trim import apply_trim_configs
import abi_sauce.streamlit_cache as streamlit_cache
from abi_sauce.streamlit_cache import (
    build_assembly_trace_row_source_for_member,
    build_parsed_batch_cache_key,
    build_selected_assembly_trace_view,
    build_trim_inputs_cache_key,
    compute_saved_assemblies_for_definitions,
    prepare_batch_for_trim_inputs,
)
from abi_sauce.trim_state import ResolvedBatchTrimInputs
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
        description="synthetic cached record",
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


def make_parsed_batch() -> ParsedBatch:
    uploads = (
        SequenceUpload(filename="left.ab1", content=b"left"),
        SequenceUpload(filename="right.ab1", content=b"right"),
        SequenceUpload(filename="short.ab1", content=b"short"),
    )
    return ParsedBatch(
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


def make_orientation_flipped_parsed_batch(
    parsed_batch: ParsedBatch,
) -> ParsedBatch:
    flipped_left_record = replace(
        parsed_batch.parsed_records["left.ab1"],
        orientation="reverse_complement",
    )
    return ParsedBatch(
        uploads=parsed_batch.uploads,
        parsed_records={
            **parsed_batch.parsed_records,
            "left.ab1": flipped_left_record,
        },
        parse_errors=dict(parsed_batch.parse_errors),
        signature=parsed_batch.signature,
    )


def test_build_trim_inputs_cache_key_normalizes_record_config_order() -> None:
    first = ResolvedBatchTrimInputs(
        default_trim_config=TrimConfig(left_trim=1),
        trim_configs_by_name={
            "b.ab1": TrimConfig(right_trim=1),
            "a.ab1": TrimConfig(left_trim=2),
        },
    )
    second = ResolvedBatchTrimInputs(
        default_trim_config=TrimConfig(left_trim=1),
        trim_configs_by_name={
            "a.ab1": TrimConfig(left_trim=2),
            "b.ab1": TrimConfig(right_trim=1),
        },
    )

    assert build_trim_inputs_cache_key(first) == build_trim_inputs_cache_key(second)


def test_build_parsed_batch_cache_key_changes_when_record_orientation_changes() -> None:
    parsed_batch = make_parsed_batch()
    flipped_parsed_batch = make_orientation_flipped_parsed_batch(parsed_batch)

    assert build_parsed_batch_cache_key(parsed_batch) != build_parsed_batch_cache_key(
        flipped_parsed_batch
    )


def test_prepare_batch_for_trim_inputs_matches_batch_service() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    resolved_trim_inputs = ResolvedBatchTrimInputs(
        default_trim_config=TrimConfig(left_trim=1),
        trim_configs_by_name={"short.ab1": TrimConfig(right_trim=2, min_length=3)},
    )

    cached_prepared_batch = prepare_batch_for_trim_inputs(
        parsed_batch,
        resolved_trim_inputs,
    )
    expected_prepared_batch = apply_trim_configs(
        parsed_batch,
        default_trim_config=resolved_trim_inputs.default_trim_config,
        trim_configs_by_name=resolved_trim_inputs.trim_configs_by_name,
    )

    assert cached_prepared_batch == expected_prepared_batch


def test_prepare_batch_for_trim_inputs_respects_orientation_changes() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    flipped_parsed_batch = make_orientation_flipped_parsed_batch(parsed_batch)
    resolved_trim_inputs = ResolvedBatchTrimInputs(default_trim_config=TrimConfig())

    original_prepared_batch = prepare_batch_for_trim_inputs(
        parsed_batch,
        resolved_trim_inputs,
    )
    flipped_prepared_batch = prepare_batch_for_trim_inputs(
        flipped_parsed_batch,
        resolved_trim_inputs,
    )

    assert (
        original_prepared_batch.trim_results["left.ab1"].record.orientation == "forward"
    )
    assert (
        flipped_prepared_batch.trim_results["left.ab1"].record.orientation
        == "reverse_complement"
    )
    assert original_prepared_batch != flipped_prepared_batch


def test_compute_saved_assemblies_for_definitions_matches_direct_service() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    prepared_batch = apply_trim_configs(parsed_batch)
    definitions = (
        AssemblyDefinition(
            assembly_id="assembly-good",
            name="Amplicon A",
            source_filenames=("left.ab1", "right.ab1"),
            config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        ),
        AssemblyDefinition(
            assembly_id="assembly-multi",
            name="Amplicon Multi",
            source_filenames=("left.ab1", "right.ab1", "short.ab1"),
            config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
            engine_kind="multi",
        ),
    )

    cached_computed_assemblies = compute_saved_assemblies_for_definitions(
        prepared_batch,
        definitions,
    )
    expected_computed_assemblies = compute_saved_assemblies(
        prepared_batch,
        definitions,
    )

    assert cached_computed_assemblies == expected_computed_assemblies
    assert tuple(cached_computed_assemblies) == ("assembly-good", "assembly-multi")


def test_compute_saved_assemblies_for_definitions_varies_with_prepared_batch_state() -> (
    None
):
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    flipped_parsed_batch = make_orientation_flipped_parsed_batch(parsed_batch)
    definitions = (
        AssemblyDefinition(
            assembly_id="assembly-good",
            name="Amplicon A",
            source_filenames=("left.ab1", "right.ab1"),
            config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        ),
    )

    original_prepared_batch = apply_trim_configs(parsed_batch)
    flipped_prepared_batch = apply_trim_configs(flipped_parsed_batch)

    original_computed_assemblies = compute_saved_assemblies_for_definitions(
        original_prepared_batch,
        definitions,
    )
    flipped_computed_assemblies = compute_saved_assemblies_for_definitions(
        flipped_prepared_batch,
        definitions,
    )

    assert (
        original_computed_assemblies["assembly-good"].result
        != flipped_computed_assemblies["assembly-good"].result
    )


def test_build_selected_assembly_trace_view_matches_pairwise_builder() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    prepared_batch = apply_trim_configs(parsed_batch)
    definition = AssemblyDefinition(
        assembly_id="assembly-good",
        name="Amplicon A",
        source_filenames=("left.ab1", "right.ab1"),
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
    )
    computed_assembly = compute_saved_assemblies_for_definitions(
        prepared_batch,
        (definition,),
    )["assembly-good"]

    cached_trace_view = build_selected_assembly_trace_view(
        prepared_batch,
        computed_assembly,
    )
    assert computed_assembly.result is not None
    assert not isinstance(computed_assembly.result, MultiAssemblyResult)

    direct_trace_view = build_pairwise_assembly_trace_view(
        result=computed_assembly.result,
        left_source_filename="left.ab1",
        left_raw_record=prepared_batch.parsed_records["left.ab1"],
        left_trim_result=prepared_batch.trim_results["left.ab1"],
        right_source_filename="right.ab1",
        right_raw_record=prepared_batch.parsed_records["right.ab1"],
        right_trim_result=prepared_batch.trim_results["right.ab1"],
    )

    assert cached_trace_view == direct_trace_view


def test_build_selected_assembly_trace_view_matches_multi_builder() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    prepared_batch = apply_trim_configs(parsed_batch)
    definition = AssemblyDefinition(
        assembly_id="assembly-multi",
        name="Amplicon Multi",
        source_filenames=("left.ab1", "right.ab1", "short.ab1"),
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
        engine_kind="multi",
    )
    computed_assembly = compute_saved_assemblies_for_definitions(
        prepared_batch,
        (definition,),
    )["assembly-multi"]

    cached_trace_view = build_selected_assembly_trace_view(
        prepared_batch,
        computed_assembly,
    )
    assert computed_assembly.result is not None
    assert isinstance(computed_assembly.result, MultiAssemblyResult)

    direct_trace_view = build_multi_assembly_trace_view(
        result=computed_assembly.result,
        raw_records_by_source_filename=prepared_batch.parsed_records,
        trim_results_by_source_filename=prepared_batch.trim_results,
    )

    assert cached_trace_view == direct_trace_view


def test_build_assembly_trace_row_source_for_member_matches_direct_builder() -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()
    prepared_batch = apply_trim_configs(parsed_batch)

    cached_row_source = build_assembly_trace_row_source_for_member(
        prepared_batch.parsed_records["right.ab1"],
        prepared_batch.trim_results["right.ab1"],
        strand="reverse_complement",
    )
    direct_row_source = build_assembly_trace_row_source(
        raw_record=prepared_batch.parsed_records["right.ab1"],
        trim_result=prepared_batch.trim_results["right.ab1"],
        strand="reverse_complement",
    )

    assert cached_row_source == direct_row_source


def test_prepare_batch_for_trim_inputs_reuses_cached_trim_results_for_unchanged_records(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    st.cache_data.clear()
    parsed_batch = make_parsed_batch()

    original_trim_sequence_record = streamlit_cache.trim_sequence_record
    trim_calls: list[tuple[str, TrimConfig]] = []

    def counting_trim_sequence_record(record: SequenceRecord, trim_config: TrimConfig):
        trim_calls.append((record.name, trim_config))
        return original_trim_sequence_record(record, trim_config)

    monkeypatch.setattr(
        streamlit_cache,
        "trim_sequence_record",
        counting_trim_sequence_record,
    )

    prepare_batch_for_trim_inputs(
        parsed_batch,
        ResolvedBatchTrimInputs(
            default_trim_config=TrimConfig(left_trim=1),
        ),
    )
    assert len(trim_calls) == 3

    trim_calls.clear()
    prepare_batch_for_trim_inputs(
        parsed_batch,
        ResolvedBatchTrimInputs(
            default_trim_config=TrimConfig(left_trim=1),
            trim_configs_by_name={"short.ab1": TrimConfig(right_trim=2)},
        ),
    )

    assert trim_calls == [("short", TrimConfig(right_trim=2))]
