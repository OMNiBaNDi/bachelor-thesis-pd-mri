"""Subject ID normalization.

The same subject shows up as S05, S5, "S 5", and S005 across different
sources. This pads trailing digits to 3 places so they all match.
"""
from __future__ import annotations

import re


_SPLIT_LETTERS_DIGITS = re.compile(r"^([A-Za-z]+)(\d+)$")


def normalize_subject_id(raw: str) -> str:
    """Pad trailing digits to 3 places: S05 -> S005, BK17 -> BK017.

    Strings that don't match <letters><digits> get whitespace stripped
    and are otherwise returned as-is. IDs already at 3+ digits pass
    through.
    """
    if not isinstance(raw, str):
        return raw
    s = raw.strip().replace(" ", "")
    m = _SPLIT_LETTERS_DIGITS.match(s)
    if not m:
        return s
    prefix, digits = m.group(1), m.group(2)
    if len(digits) < 3:
        digits = digits.zfill(3)
    return prefix + digits
