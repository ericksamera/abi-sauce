from __future__ import annotations

from datetime import UTC, datetime


def utcnow() -> datetime:
    """Timezone-aware 'now' in UTC."""
    return datetime.now(UTC)
