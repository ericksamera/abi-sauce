from __future__ import annotations

from uuid import uuid4

from abi_sauce.models import AssetBase, Sample
from abi_sauce.services.file_manager import FileManager


class SampleManager:
    """Owns user-visible 'samples' built on top of raw assets (non-destructive)."""

    def __init__(self, file_manager: FileManager) -> None:
        self.fm = file_manager
        self._samples: dict[str, Sample] = {}

    # --- create/import ---
    def import_bytes(self, *, name: str, raw: bytes) -> list[str]:
        """Parse bytes via FileManager, create one Sample per new asset."""
        new_asset_ids = self.fm.add_bytes(name=name, raw=raw)
        sample_ids: list[str] = []
        for aid in new_asset_ids:
            a = self.fm.get(aid)
            sid = uuid4().hex
            sample = Sample(
                id=sid,
                name=a.name,
                asset_ids=[aid],
                primary_asset_id=aid,
            )
            self._samples[sid] = sample
            sample_ids.append(sid)
        return sample_ids

    # --- queries ---
    def get(self, sample_id: str) -> Sample:
        return self._samples[sample_id]

    def list(self) -> list[Sample]:
        return sorted(self._samples.values(), key=lambda s: s.created_at, reverse=True)

    def assets_view(self) -> dict[str, AssetBase]:
        # convenience map for render/export
        return {
            aid: self.fm.get(aid)
            for aid in {aid for s in self._samples.values() for aid in s.asset_ids}
        }

    # --- edits ---
    def rename(self, sample_id: str, new_name: str) -> None:
        self._samples[sample_id].name = new_name

    def set_primary(self, sample_id: str, asset_id: str) -> None:
        s = self._samples[sample_id]
        if asset_id in s.asset_ids:
            s.primary_asset_id = asset_id

    def set_sequence_override(self, sample_id: str, seq: str | None) -> None:
        """Set or clear the editable sequence (None clears; empty string also clears)."""
        if seq is None:
            self._samples[sample_id].sequence_override = None
            return
        sanitized = seq.replace("\n", "").replace("\r", "")
        self._samples[sample_id].sequence_override = sanitized or None

    def set_feature_overrides(self, sample_id: str, feats: list[dict] | None) -> None:
        """None => inherit features from primary asset; list => explicit override."""
        self._samples[sample_id].feature_overrides = (
            None if feats is None else list(feats)
        )

    def set_tags(self, sample_id: str, tags: list[str] | None) -> None:
        """Replace tags (whitespace-trimmed, empties removed)."""
        self._samples[sample_id].tags = [t.strip() for t in (tags or []) if t.strip()]

    # --- export ---
    def fasta(self, sample_id: str) -> str | None:
        s = self._samples[sample_id]
        return s.to_fasta(self.fm._assets)
