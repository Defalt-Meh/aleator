"""Aleator: Monte Carlo, molecular dynamics, and porous-material geometry
analysis on a shared periodic-system core.

This build exposes scaffolding only — no physics is implemented yet. See the
project README for current status.
"""

from __future__ import annotations

from ._aleator import __version__ as __version__
from ._aleator import vector_sum as vector_sum

__all__ = ["__version__", "vector_sum"]
