from __future__ import annotations

from hashlib import md5


def md5_hex(data: bytes) -> str:
    """Hex MD5 for consistent, typed hashing."""
    return md5(data).hexdigest()
