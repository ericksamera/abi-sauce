from __future__ import annotations


class AbiSauceError(Exception):
    """Base exception for abi_sauce."""


class ParseError(AbiSauceError):
    """Base exception for sequence parsing failures."""


class AbiParseError(ParseError):
    """Raised when an ABI file cannot be parsed."""
