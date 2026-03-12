from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class SequenceUpload:
    """Framework-agnostic uploaded sequence file."""

    filename: str
    content: bytes

    @property
    def suffix(self) -> str:
        """Return the lowercase file extension without the leading dot."""
        return Path(self.filename).suffix.lower().lstrip(".")

    @property
    def size_bytes(self) -> int:
        """Return the upload size in bytes."""
        return len(self.content)


@dataclass(slots=True)
class TraceData:
    """Chromatogram-level information for trace-based formats such as ABI."""

    channels: dict[str, list[int]]
    base_positions: list[int] = field(default_factory=list)
    channel_order: str | None = None


@dataclass(slots=True)
class SequenceRecord:
    """Normalized application-level sequence record."""

    record_id: str
    name: str
    description: str
    sequence: str
    source_format: str
    qualities: list[int] | None = None
    trace_data: TraceData | None = None
    annotations: dict[str, Any] = field(default_factory=dict)
