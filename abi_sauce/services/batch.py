from __future__ import annotations

from collections.abc import Iterable

from abi_sauce.exceptions import AbiParseError
from abi_sauce.models import SequenceUpload
from abi_sauce.parsers.abi import parse_ab1_upload
from abi_sauce.services.batch_export import (
    BatchDownloadArtifact,
    BatchExportSelection,
    prepare_batch_download,
    select_batch_export,
)
from abi_sauce.services.batch_parse import (
    BatchSignature,
    ParsedBatch,
    UploadedFileLike,
    build_batch_signature,
    normalize_uploaded_files,
    parse_uploaded_batch as _parse_uploaded_batch,
    parse_uploads as _parse_uploads,
    replace_parsed_batch_record,
)
from abi_sauce.services.batch_trim import (
    PreparedBatch,
    apply_trim_config,
    apply_trim_configs,
)


def parse_uploaded_batch(uploaded_files: Iterable[UploadedFileLike]) -> ParsedBatch:
    """Normalize and parse one batch of uploaded ABI files."""
    return _parse_uploaded_batch(
        uploaded_files,
        parse_upload=parse_ab1_upload,
    )


def parse_uploads(uploads: Iterable[SequenceUpload]) -> ParsedBatch:
    """Parse normalized uploads into records plus parse failures."""
    return _parse_uploads(
        uploads,
        parse_upload=parse_ab1_upload,
    )


__all__ = [
    "AbiParseError",
    "BatchDownloadArtifact",
    "BatchExportSelection",
    "BatchSignature",
    "ParsedBatch",
    "PreparedBatch",
    "UploadedFileLike",
    "apply_trim_config",
    "apply_trim_configs",
    "build_batch_signature",
    "normalize_uploaded_files",
    "parse_ab1_upload",
    "parse_uploaded_batch",
    "parse_uploads",
    "prepare_batch_download",
    "replace_parsed_batch_record",
    "select_batch_export",
]
