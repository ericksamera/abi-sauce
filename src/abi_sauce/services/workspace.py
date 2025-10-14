from __future__ import annotations

import contextlib
import json
from datetime import datetime

from abi_sauce.models import Sample, SequenceAsset, TraceAsset
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager


class WorkspaceManager:
    """Serialize/restore assets + samples for a running session.
    This touches FileManager internals by design (UI decides if/when to persist)."""

    def export_json(self, fm: FileManager, sm: SampleManager, indent: int = 2) -> str:
        # assets -> plain dicts
        assets = []
        for a in fm._assets.values():
            if isinstance(a, SequenceAsset):
                payload = {
                    "type": "sequence",
                    "data": {
                        "id": a.id,
                        "name": a.name,
                        "ext": a.ext,
                        "size": a.size,
                        "checksum_md5": a.checksum_md5,
                        "sequence": a.sequence,
                        "description": a.description,
                        "length": a.length,
                        "meta": a.meta,
                        "features": a.features,
                        "created_at": a.created_at.isoformat(),
                    },
                }
            elif isinstance(a, TraceAsset):
                payload = {
                    "type": "trace",
                    "data": {
                        "id": a.id,
                        "name": a.name,
                        "ext": a.ext,
                        "size": a.size,
                        "checksum_md5": a.checksum_md5,
                        "sequence": a.sequence,
                        "qualities": a.qualities,
                        "base_positions": a.base_positions,
                        "channels": a.channels,
                        "meta": a.meta,
                        "created_at": a.created_at.isoformat(),
                    },
                }
            else:
                continue
            assets.append(payload)

        # samples -> dicts
        samples = []
        for s in sm._samples.values():
            samples.append(
                {
                    "id": s.id,
                    "name": s.name,
                    "asset_ids": list(s.asset_ids),
                    "description": s.description,
                    "tags": list(s.tags),
                    "primary_asset_id": s.primary_asset_id,
                    "sequence_override": s.sequence_override,
                    "feature_overrides": s.feature_overrides,
                    "created_at": s.created_at.isoformat(),
                }
            )

        return json.dumps({"assets": assets, "samples": samples}, indent=indent)

    def import_json(
        self,
        fm: FileManager,
        sm: SampleManager,
        payload: str,
        *,
        clear_first: bool = True,
    ) -> tuple[int, int]:
        """Clear managers and restore a snapshot. Returns (n_assets, n_samples)."""
        snapshot = json.loads(payload)
        if clear_first:
            fm.clear()
            sm._samples.clear()

        # restore assets
        by_checksum: dict[str, list] = {}
        for item in snapshot.get("assets", []):
            typ = item["type"]
            d = item["data"]
            if typ == "sequence":
                a = SequenceAsset(
                    id=d["id"],
                    name=d["name"],
                    ext=d["ext"],
                    size=d["size"],
                    checksum_md5=d["checksum_md5"],
                    sequence=d["sequence"],
                    description=d.get("description", ""),
                    length=d.get("length", 0),
                    meta=d.get("meta", {}),
                    features=d.get("features", []),
                )
            elif typ == "trace":
                a = TraceAsset(
                    id=d["id"],
                    name=d["name"],
                    ext=d["ext"],
                    size=d["size"],
                    checksum_md5=d["checksum_md5"],
                    sequence=d.get("sequence"),
                    qualities=d.get("qualities"),
                    base_positions=d.get("base_positions"),
                    channels=d.get("channels", {}),
                    meta=d.get("meta", {}),
                )
            else:
                continue
            # created_at is init=False, set post-init
            with contextlib.suppress(Exception):
                a.created_at = datetime.fromisoformat(d["created_at"])
            fm._assets[a.id] = a
            by_checksum.setdefault(a.checksum_md5, []).append(a.id)
        fm._by_checksum = by_checksum

        # restore samples
        for s in snapshot.get("samples", []):
            sample = Sample(
                id=s["id"],
                name=s["name"],
                asset_ids=list(s.get("asset_ids") or []),
                description=s.get("description", ""),
                tags=list(s.get("tags") or []),
                primary_asset_id=s.get("primary_asset_id"),
            )
            sample.sequence_override = s.get("sequence_override")
            sample.feature_overrides = s.get("feature_overrides")
            with contextlib.suppress(Exception):
                sample.created_at = datetime.fromisoformat(s["created_at"])
            sm._samples[sample.id] = sample

        return len(fm._assets), len(sm._samples)
