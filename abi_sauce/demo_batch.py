from __future__ import annotations

from pathlib import Path

from abi_sauce.models import SequenceUpload
from abi_sauce.services.batch import parse_uploads
from abi_sauce.upload_state import set_active_parsed_batch
from abi_sauce.viewer_state import clear_viewer_session_state
from abi_sauce.assembly_state import clear_assembly_session_state


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def demo_sample_path() -> Path:
    return _repo_root() / "tests" / "fixtures" / "example.ab1"


def can_load_demo_sample() -> bool:
    return demo_sample_path().exists()


def load_demo_sample(session_state) -> None:
    sample_path = demo_sample_path()
    upload = SequenceUpload(
        filename=sample_path.name,
        content=sample_path.read_bytes(),
    )
    clear_viewer_session_state(session_state)
    clear_assembly_session_state(session_state)
    set_active_parsed_batch(session_state, parse_uploads((upload,)))
