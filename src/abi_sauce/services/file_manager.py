from __future__ import annotations

from abi_sauce.io import readers
from abi_sauce.io.detect import identify
from abi_sauce.models import AssetBase, AssetKind
from abi_sauce.services.errors import UnsupportedFileType
from abi_sauce.utils.hashing import md5_hex


class FileManager:
    """In-memory registry of uploaded assets (sequence + trace).

    Dedup policy: checksum-only (identical bytes -> skip).
    """

    def __init__(self) -> None:
        self._assets: dict[str, AssetBase] = {}
        self._by_checksum: dict[str, list[str]] = {}

    def add_bytes(self, *, name: str, raw: bytes) -> list[str]:
        size = len(raw)
        checksum = md5_hex(raw)
        kind, norm_ext = identify(name, raw)

        if kind == "unknown":
            raise UnsupportedFileType(f"Unsupported or unknown file type for: {name}")

        # de-dup check
        if self._by_checksum.get(checksum):
            return []  # nothing added

        reader = readers.READERS.get(kind)
        if not reader:
            raise UnsupportedFileType(f"No reader registered for: {kind}")

        # genbank/ape need an extension hint; others just ignore the kwarg
        assets = reader(name, raw, size, checksum, ext_hint=norm_ext)

        new_ids: list[str] = []
        for asset in assets:
            self._assets[asset.id] = asset
            self._by_checksum.setdefault(checksum, []).append(asset.id)
            new_ids.append(asset.id)
        return new_ids

    # ---- queries ----
    def get(self, asset_id: str) -> AssetBase:
        return self._assets[asset_id]

    def list(self, kind: AssetKind | None = None) -> list[AssetBase]:
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
