"""Tests for kaitiaki.utils.transforms: small axis-aware transform
functions used by plot()/kippenhahn_diagram() to optionally log/unlog or
take the absolute value of one axis at a time.
"""
import numpy as np
import pytest

from kaitiaki.utils import transforms


@pytest.mark.parametrize("axis,expected", [("x", [1, 2, 3]),
                                           ("y", [10, 100, 1000])])
def test_unlog_y_only_transforms_y(axis, expected):
    arr = np.array([1, 2, 3])

    result = transforms.unlog_y(axis, arr)

    np.testing.assert_allclose(result, expected)


@pytest.mark.parametrize("axis,expected", [("y", [1, 2, 3]),
                                           ("x", [10, 100, 1000])])
def test_unlog_x_only_transforms_x(axis, expected):
    arr = np.array([1, 2, 3])

    result = transforms.unlog_x(axis, arr)

    np.testing.assert_allclose(result, expected)


def test_log10_y_and_x():
    arr = np.array([1, 10, 100])

    np.testing.assert_allclose(transforms.log10_y("y", arr), [0, 1, 2])
    np.testing.assert_allclose(transforms.log10_y("x", arr), arr)
    np.testing.assert_allclose(transforms.log10_x("x", arr), [0, 1, 2])
    np.testing.assert_allclose(transforms.log10_x("y", arr), arr)


def test_loge_y_and_x():
    arr = np.array([1, np.e, np.e**2])

    np.testing.assert_allclose(transforms.loge_y("y", arr), [0, 1, 2])
    np.testing.assert_allclose(transforms.loge_y("x", arr), arr)
    np.testing.assert_allclose(transforms.loge_x("x", arr), [0, 1, 2])
    np.testing.assert_allclose(transforms.loge_x("y", arr), arr)


def test_unlog_both_ignores_axis():
    arr = np.array([1, 2, 3])

    np.testing.assert_allclose(transforms.unlog_both("x", arr), 10.0**arr)
    np.testing.assert_allclose(transforms.unlog_both("y", arr), 10.0**arr)


def test_absval_only_transforms_y():
    arr = np.array([-1, -2, 3])

    np.testing.assert_allclose(transforms.absval("y", arr), [1, 2, 3])
    np.testing.assert_allclose(transforms.absval("x", arr), arr)


def test_null_is_identity():
    arr = np.array([-1, 2, -3])

    np.testing.assert_allclose(transforms.null("x", arr), arr)
    np.testing.assert_allclose(transforms.null("y", arr), arr)
