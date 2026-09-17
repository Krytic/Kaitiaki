"""Tests for kaitiaki.kipp.decode: turning a plot-file row's 12 boundary
values into mass-coordinate Interval objects.

See the module's own docstring for the physics (rad<->semi<->conv state
walk). These tests construct small, hand-derived boundary sequences rather
than real STARS output, so the expected states can be verified by tracing
the state machine by hand.
"""
import pytest

from kaitiaki.kipp.decode import (
    Interval,
    DecodeError,
    decode_row,
    decode_all,
    is_truncated,
    unknown_regions,
)


def test_interval_is_a_frozen_dataclass():
    iv = Interval(0.0, 5.0, "conv")

    assert (iv.lo, iv.hi, iv.kind) == (0.0, 5.0, "conv")
    with pytest.raises(Exception):
        iv.lo = 1.0  # frozen -> should refuse mutation


def test_decode_row_all_padding_is_purely_radiative():
    values = [0.0] * 12  # every slot is padding -> no boundaries at all

    intervals = decode_row(values, m_total=10.0)

    assert intervals == []


def test_decode_row_single_boundary_gives_conv_then_semi():
    # One conv/semi boundary (negative) at mass coordinate 3, in a 10
    # Msun star. Starting from 'conv' (the walk's tie-break preference),
    # this must produce a conv interval [0, 3) and a semi interval
    # [3, 10).
    values = [-3.0] + [0.0] * 11

    intervals = decode_row(values, m_total=10.0)

    assert intervals == [
        Interval(0.0, 3.0, "conv"),
        Interval(3.0, 10.0, "semi"),
    ]


def test_decode_row_raises_decodeerror_for_illegal_sequence():
    # +, -, + at distinct (non-tied) mass coordinates: every one of the
    # three possible start states (rad/conv/semi) hits an illegal
    # transition somewhere in this specific sequence (traced by hand in
    # the module's algorithm: rad and conv both fail at the first or
    # third boundary, and semi fails at the second).
    values = [1.0, -2.0, 3.0] + [0.0] * 9

    with pytest.raises(DecodeError):
        decode_row(values, m_total=10.0)


def test_is_truncated_true_when_all_twelve_slots_used():
    values = list(range(1, 13))  # 12 distinct non-padding magnitudes

    assert is_truncated(values, m_total=20.0) is True


def test_is_truncated_false_with_padding_present():
    values = [1.0, 2.0] + [0.0] * 10

    assert is_truncated(values, m_total=20.0) is False


def test_decode_all_skip_records_bad_rows_without_raising():
    good_row = [-3.0] + [0.0] * 11
    bad_row = [1.0, -2.0, 3.0] + [0.0] * 9

    results, bad_rows = decode_all([good_row, bad_row], M=[10.0, 10.0],
                                   on_error="skip")

    assert bad_rows == [1]
    assert results[0] != []
    assert results[1] == []


def test_decode_all_raise_propagates_decodeerror_with_row_index():
    bad_row = [1.0, -2.0, 3.0] + [0.0] * 9

    with pytest.raises(DecodeError, match="row 0"):
        decode_all([bad_row], M=[10.0], on_error="raise")


def test_decode_all_rejects_bad_on_error_value():
    with pytest.raises(ValueError):
        decode_all([], M=[], on_error="explode")


def test_decode_all_restores_truncated_envelope_from_conv_env():
    # All 12 slots used -> truncated=True, so decode_all should append a
    # conv interval reaching from wherever the boundaries stop out to
    # conv_env (below the surface), using the pattern [+,+,-,-] repeated
    # (always walkable starting from the 'semi' state -- see decode.py's
    # docstring for the state machine rules).
    signs = [1, 1, -1, -1] * 3
    magnitudes = range(1, 13)
    row = [sign * mag for sign, mag in zip(signs, magnitudes)]
    m_total = 20.0
    conv_env = 15.0

    results, bad_rows = decode_all([row], M=[m_total], conv_env=[conv_env])

    assert bad_rows == []
    intervals = results[0]
    assert intervals  # something was decoded
    restored = intervals[-1]
    assert restored.kind == "conv"
    assert restored.hi == pytest.approx(m_total)


def test_unknown_regions_normal_row_is_none():
    row = [-3.0] + [0.0] * 11

    regions = unknown_regions([row], M=[10.0])

    assert regions == [None]


def test_unknown_regions_truncated_row_reports_range():
    signs = [1, 1, -1, -1] * 3
    magnitudes = range(1, 13)
    row = [sign * mag for sign, mag in zip(signs, magnitudes)]
    m_total = 20.0

    regions = unknown_regions([row], M=[m_total], conv_env=[15.0])

    assert regions[0] is not None
    lo, hi = regions[0]
    assert lo == pytest.approx(12.0)  # outermost reported boundary
    assert hi == pytest.approx(15.0)  # clipped to conv_env
