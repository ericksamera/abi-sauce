from __future__ import annotations

from pathlib import Path

import pytest

from abi_sauce.models import SequenceUpload

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def real_ab1_path() -> Path:
    path = FIXTURES_DIR / "example.ab1"
    if not path.exists():
        pytest.skip("Missing test fixture: tests/fixtures/example.ab1")
    return path


@pytest.fixture(scope="session")
def real_ab1_upload(real_ab1_path: Path) -> SequenceUpload:
    return SequenceUpload(
        filename=real_ab1_path.name,
        content=real_ab1_path.read_bytes(),
    )
