"""Smoke tests for the Python packaging layer. No physics is tested here —
only that the extension module builds, imports, and the one trivial bound
function behaves as documented."""

import aleator


def test_version_is_a_nonempty_string():
    assert isinstance(aleator.__version__, str)
    assert aleator.__version__


def test_vector_sum_matches_python_sum():
    values = [1, 2, 3, 45, -7]
    assert aleator.vector_sum(values) == sum(values)


def test_vector_sum_empty():
    assert aleator.vector_sum([]) == 0
