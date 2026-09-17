"""Tests for kaitiaki.kipp.rasterise: time_edges and rasterise.

Pure numpy functions with no global state -- see the module docstring.
"""
import numpy as np
import pytest

from kaitiaki.kipp.decode import Interval
from kaitiaki.kipp.rasterise import time_edges, rasterise


def test_time_edges_basic_increasing_sequence():
    edges = time_edges([0.0, 1.0, 2.0, 4.0])

    assert len(edges) == 5
    # every original point should lie strictly between its two edges
    assert edges[0] < 0.0 < edges[1]
    assert edges[1] < 1.0 < edges[2]
    assert edges[2] < 2.0 < edges[3]
    assert edges[3] < 4.0 < edges[4]


def test_time_edges_single_point():
    edges = time_edges([5.0])

    np.testing.assert_allclose(edges, [4.5, 5.5])


def test_time_edges_strictly_decreasing_sequence_is_allowed():
    edges = time_edges([3.0, 2.0, 1.0])

    assert len(edges) == 4
    assert edges[0] > 3.0
    assert edges[-1] < 1.0


def test_time_edges_rejects_empty_input():
    with pytest.raises(ValueError):
        time_edges([])


def test_time_edges_rejects_non_monotone_input():
    with pytest.raises(ValueError):
        time_edges([0.0, 1.0, 0.5, 2.0])


def test_time_edges_rejects_constant_input():
    with pytest.raises(ValueError):
        time_edges([1.0, 1.0, 1.0])


def test_rasterise_marks_cells_inside_interval():
    intervals_per_model = [[Interval(2.0, 5.0, "conv")]]

    m_edges, conv, semi = rasterise(intervals_per_model, m_max=10.0,
                                    n_mass=10)

    assert m_edges[0] == 0.0
    assert m_edges[-1] == 10.0
    assert conv.shape == (1, 10)
    assert semi.shape == (1, 10)
    # cell centres are at 0.5, 1.5, ..., 9.5 -> [2, 5) covers centres
    # 2.5, 3.5, 4.5, i.e. columns 2, 3, 4
    assert conv[0].tolist() == [False, False, True, True, True,
                                False, False, False, False, False]
    assert not semi.any()


def test_rasterise_semi_kind_goes_to_semi_grid():
    intervals_per_model = [[Interval(0.0, 3.0, "semi")]]

    _, conv, semi = rasterise(intervals_per_model, m_max=10.0, n_mass=10)

    assert not conv.any()
    assert semi[0, :3].all()
    assert not semi[0, 3:].any()


def test_rasterise_skips_degenerate_interval():
    intervals_per_model = [[Interval(5.0, 5.0, "conv")]]  # hi <= lo

    _, conv, semi = rasterise(intervals_per_model, m_max=10.0, n_mass=10)

    assert not conv.any()
    assert not semi.any()


def test_rasterise_handles_multiple_models_independently():
    intervals_per_model = [
        [Interval(0.0, 2.0, "conv")],
        [],
        [Interval(8.0, 10.0, "semi")],
    ]

    _, conv, semi = rasterise(intervals_per_model, m_max=10.0, n_mass=10)

    assert conv[0].any()
    assert not conv[1].any() and not semi[1].any()
    assert semi[2].any()
