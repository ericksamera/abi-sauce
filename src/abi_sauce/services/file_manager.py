from __future__ import annotations
from typing import Dict, List, Optional
import hashlib

from abi_sauce.models import AssetBase, AssetKind
from abi_sauce.io.detect import identify
from abi_sauce.io import readers


class FileManager:
    """In-memory registry of uploaded assets (sequence + trace).

    - De-duplicates on checksum+name.
    - Splits multi-record inputs into multiple assets.
    - Leaves persistence policy to the UI (e.g., session or disk).
    """

    def __init__(self) -> None:
        self._assets: Dict[str, AssetBase] = {}
        self._by_checksum: Dict[str, List[str]] = {}

    @staticmethod
    def _checksum_md5(data: bytes) -> str:
        return hashlib.md5(data).hexdigest()

    def add_bytes(self, *, name: str, raw: bytes) -> List[str]:
        size = len(raw)
        checksum = self._checksum_md5(raw)
        kind, norm_ext = identify(name, raw)

        # de-dup check (exact same filename+checksum)
        existing_ids = self._by_checksum.get(checksum, [])
        if existing_ids:
            return []  # nothing added

        new_ids: List[str] = []
        if kind == "fasta":
            assets = readers.read_fasta(name, raw, size, checksum)
        elif kind in ("genbank", "ape"):
            assets = readers.read_genbank_like(
                name, raw, size, checksum, ext_hint=norm_ext
            )
        elif kind == "ab1":
            assets = readers.read_ab1(name, raw, size, checksum)
        else:
            raise ValueError(f"Unsupported or unknown file type for: {name}")

        for asset in assets:
            self._assets[asset.id] = asset
            self._by_checksum.setdefault(checksum, []).append(asset.id)
            new_ids.append(asset.id)
        return new_ids

    # ---- queries ----
    def get(self, asset_id: str) -> AssetBase:
        return self._assets[asset_id]

    def list(self, kind: Optional[AssetKind] = None) -> List[AssetBase]:
        items = list(self._assets.values())
        if kind:
            items = [a for a in items if a.kind == kind]
        # newest first
        return sorted(items, key=lambda a: a.created_at, reverse=True)

    def remove(self, asset_id: str) -> None:
        asset = self._assets.pop(asset_id, None)
        if not asset:
            return
        # cleanup checksum index
        ids = self._by_checksum.get(asset.checksum_md5, [])
        if asset_id in ids:
            ids.remove(asset_id)
        if not ids:
            self._by_checksum.pop(asset.checksum_md5, None)

    def clear(self) -> None:
        self._assets.clear()
        self._by_checksum.clear()
