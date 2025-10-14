# src/abi_sauce/__init__.py
from __future__ import annotations

__all__ = ["__version__"]
__version__ = "0.2.0"

# Re-export commonly used types (nice for downstream users)
from .types import (  # noqa: F401
    BaseChar,
    BaseIndex,
    Channel,
    ChannelsMap,
    PerBase,
    PerBaseMap,
)
