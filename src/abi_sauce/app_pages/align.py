# Thin wrapper to keep the public entrypoint stable.
from __future__ import annotations

from abi_sauce.ui.alignment_view import align_page

# Re-export so tooling doesn't strip it as unused.
__all__ = ["align_page"]
