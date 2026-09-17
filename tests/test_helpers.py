"""Tests for kaitiaki.helpers: the Range class, find_runs, and the
Roche-lobe / contact-phase geometry helpers.
"""
import numpy as np
import pytest

from kaitiaki.helpers import RL, roche_lobes, contact_phases, find_runs, Range


def test_RL_equal_mass_matches_known_eggleton_value():
    # For q = m2/m1 = 1, Eggleton (1983)'s approximation gives R_L1/a ~ 0.38.
    assert RL(1.0, 1.0, lobe="L1") == pytest.approx(0.3789, abs=1e-3)


def test_RL_L2_is_larger_than_L1():
    q = 0.5
    a = 3.0

    assert RL(q, a, lobe="L2") > RL(q, a, lobe="L1")


def test_RL_rejects_unknown_lobe():
    with pytest.raises(AssertionError):
        RL(1.0, 1.0, lobe="L3")


class _FakeStar:
    """Minimal stand-in for a kaitiaki.file.plot object.

    Only implements .get(key), which is all that
    kaitiaki.helpers.roche_lobes and kaitiaki.helpers.contact_phases
    actually use.
    """
    def __init__(self, **columns):
        self._columns = {k: np.asarray(v) for k, v in columns.items()}

    def get(self, key):
        return self._columns[key]


def test_roche_lobes_returns_L1_and_L2():
    donor = _FakeStar(M=np.array([1.0]), **{"log(R)": np.array([0.0])})
    accretor = _FakeStar(M=np.array([1.0]), a=np.array([10.0]))

    lobes = roche_lobes(donor, accretor)

    assert set(lobes.keys()) == {"L1", "L2"}
    assert lobes["L1"] > 0
    assert lobes["L2"] > lobes["L1"]


def test_contact_phases_detects_overflow():
    # Both stars are given a huge radius (10**2 = 100 Rsun) at a tiny
    # orbital separation, guaranteeing both fill (and overflow) their
    # Roche lobes at every timestep.
    donor = _FakeStar(M=np.array([1.0, 1.0]),
                      **{"log(R)": np.array([2.0, 2.0]),
                         "a": np.array([0.1, 0.1])})
    accretor = _FakeStar(M=np.array([1.0, 1.0]),
                         **{"log(R)": np.array([2.0, 2.0]),
                            "a": np.array([0.1, 0.1])})

    in_contact = contact_phases(donor, accretor)

    assert in_contact.all()


def test_find_runs_basic():
    values, starts, lengths = find_runs(np.array([1, 1, 2, 2, 2, 3]))

    assert values.tolist() == [1, 2, 3]
    assert starts.tolist() == [0, 2, 5]
    assert lengths.tolist() == [2, 3, 1]


def test_find_runs_empty_array():
    values, starts, lengths = find_runs(np.array([]))

    assert len(values) == 0
    assert len(starts) == 0
    assert len(lengths) == 0


def test_find_runs_rejects_2d_array():
    with pytest.raises(ValueError):
        find_runs(np.zeros((2, 2)))


@pytest.mark.parametrize("inclusivity,lo_in,hi_in", [
    ("[]", True, True),
    ("()", False, False),
    ("(]", False, True),
    ("[)", True, False),
])
def test_range_inclusivity(inclusivity, lo_in, hi_in):
    r = Range(0, 10, inclusivity=inclusivity)

    assert (0 in r) is lo_in
    assert (10 in r) is hi_in
    assert 5 in r
    assert -1 not in r
    assert 11 not in r


def test_range_getitem():
    r = Range(1, 9)

    assert r[0] == 1
    assert r[1] == 9
    with pytest.raises(ValueError):
        r[2]


def test_range_setitem_requires_mutable():
    r = Range(1, 9, mutable=False)

    with pytest.raises(TypeError):
        r[0] = 5

    mutable_r = Range(1, 9, mutable=True)
    mutable_r[0] = 5

    assert mutable_r[0] == 5


def test_range_str():
    assert str(Range(1, 9, inclusivity="[)")) == "[1, 9)"
