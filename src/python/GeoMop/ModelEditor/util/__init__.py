"""
Utilities package for ModelEditor.

This package can only depend on ``common_util`` to avoid circular
references between modules.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

from .enums import CursorType, KeyType, PosType
from .util import TextValue
from .locators import Span, Position
