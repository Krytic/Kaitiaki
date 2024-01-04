import numpy as np


def unlog_y(ax, arr):
    return 10**arr if ax == 'y' else arr


def unlog_x(ax, arr):
    return 10**arr if ax == 'x' else arr


def log10_y(ax, arr):
    return np.log10(arr) if ax == 'y' else arr


def log10_x(ax, arr):
    return np.log10(arr) if ax == 'x' else arr


def loge_y(ax, arr):
    return np.log(arr) if ax == 'y' else arr


def loge_x(ax, arr):
    return np.log(arr) if ax == 'x' else arr


def unlog_both(ax, arr):
    return 10**arr


def absval(ax, arr):
    return np.abs(arr) if ax == 'y' else arr


def null(ax, arr):
    return arr