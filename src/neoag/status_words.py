"""Shared polarity helpers for free-text status fields.

Audit references (B-E 136-finding review): B5-01, B5-02, D1-03, D1-04 and
review principle F3 ("exact enums and NOT-prefix-first handling ... audit
every ``in status`` expression, especially hard-exclusion branches").

The hazard this module removes: code tested for support with substring
membership such as ``"PASS" in status`` or ``"ALT_SUPPORTED" in status``.
A negated token contains its positive token -- ``NOT_PASSED`` contains
``PASS``, ``NOT_ALT_SUPPORTED`` contains ``ALT_SUPPORTED``, ``UNPHASED``
contains ``PHASED``, ``NOT_SOMATIC`` contains ``SOMATIC`` -- so an explicit
negative was previously read as support, in hard-exclusion branches as well
as support branches.

Every positive-token test goes through :func:`contains_positive`: an
occurrence of the token counts only when it is not immediately preceded by
a negation prefix (NOT_, NO_, NON_, UN).  Whole-value tests use
:func:`starts_negated`, which returns False for negations that mean
"evidence absent" (NOT_AVAILABLE, NOT_APPLICABLE, ...) so those keep their
unassessed semantics instead of becoming a negative verdict.
"""

from __future__ import annotations

NEGATION_PREFIXES = ("NOT_", "NO_", "NON_", "UN")

# Whole-value negations that mean "this evidence does not exist", not
# "this evidence is negative"; they must not flip a verdict either way.
MISSING_NEGATIONS = frozenset(
    {
        "NOT_AVAILABLE",
        "NOT_APPLICABLE",
        "NOT_PROVIDED",
        "NOT_SUPPLIED",
        "NOT_RUN",
        "NOT_CALLED",
        "NOT_DONE",
    }
)


def normalize_status(value: object) -> str:
    """Upper-case, collapse whitespace to single underscores."""

    text = " ".join(str(value if value is not None else "").strip().upper().split())
    return text.replace(" ", "_")


def contains_positive(text: object, token: object) -> bool:
    """True when ``token`` occurs in ``text`` in a non-negated position.

    An occurrence immediately preceded by a negation prefix (NOT_, NO_,
    NON_, UN) does not count: "NOT_ALT_SUPPORTED" is not ALT_SUPPORTED,
    "UNSUPPORTED" is not SUPPORTED, "UNPHASED" is not PHASED.
    """

    up = normalize_status(text)
    tok = normalize_status(token)
    if not up or not tok:
        return False
    start = 0
    while True:
        index = up.find(tok, start)
        if index < 0:
            return False
        if not any(up[:index].endswith(prefix) for prefix in NEGATION_PREFIXES):
            return True
        start = index + len(tok)


def starts_negated(text: object) -> bool:
    """True when the whole value is an explicit negation (e.g. NOT_PASSED).

    Values that merely state absence (NOT_AVAILABLE, NOT_APPLICABLE, ...)
    are not treated as negations so callers can keep unassessed semantics.
    """

    up = normalize_status(text)
    if not up or up in MISSING_NEGATIONS:
        return False
    return any(up.startswith(prefix) for prefix in NEGATION_PREFIXES)
