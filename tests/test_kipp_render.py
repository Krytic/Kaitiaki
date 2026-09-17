"""Tests for kaitiaki.kipp.render.plot_kippenhahn.

Uses the Agg backend (forced in conftest.py) so nothing tries to open a
window, and relies on conftest.py's autouse `_close_all_figures` fixture
to clean up every figure a test creates.
"""
import matplotlib.pyplot as plt
import numpy as np
import pytest

from kaitiaki.kipp.decode import Interval
from kaitiaki.kipp.render import plot_kippenhahn


def _toy_data(n=6):
    return {
        "timestep": np.arange(n, dtype=float),
        "age": np.linspace(0.0, 100.0, n),
        "M": np.linspace(20.0, 19.0, n),
        "He_core": np.linspace(0.0, 2.0, n),
        "CO_core": np.linspace(0.0, 1.0, n),
        "conv": np.zeros((n, 12)),
        "conv_env": np.full(n, np.nan),
    }


def test_plot_kippenhahn_invalid_xaxis_raises():
    with pytest.raises(ValueError):
        plot_kippenhahn(_toy_data(), xaxis="bogus")