"""Tests for kaitiaki.classifier -- the parts that don't require a full
STARS `out`/`plot` file pair.

Notes:
    classifier.go() (the main model-classification entry point) reads a
    real out/plot file pair via kaitiaki.file.out / kaitiaki.file.plot and
    is a good candidate for a follow-up integration-style test built on a
    real out/plot fixture pair; that's out of scope for this initial pass,
    which covers the pure/self-contained pieces of the module.
"""
import numpy as np
import pytest

import kaitiaki


def test_to_str_known_codes():
    assert kaitiaki.classify.to_str(1.1) == "SNe"
    assert kaitiaki.classify.to_str(7.0) == "Too old"
    assert kaitiaki.classify.to_str(0.0) == "Unclassifiable"


def test_to_str_unknown_code_raises_keyerror():
    with pytest.raises(KeyError):
        kaitiaki.classify.to_str(999.9)


def test_strings_is_a_dict_of_known_outcomes():
    strings = kaitiaki.classify.strings()

    assert isinstance(strings, dict)
    assert strings[1.1] == "SNe"
    assert -2.0 in strings


def test_compute_separation_matches_keplers_third_law():
    # A 1+1 Msun binary with a 1-day period has a separation of a few
    # solar radii; this is mainly a smoke test that the vectorised
    # implementation runs and returns a sane order of magnitude.
    a = kaitiaki.classify.compute_separation(1.0, 1.0, 1.0)

    assert 3.0 < float(a) < 6.0


def test_compute_separation_is_vectorized():
    a = kaitiaki.classify.compute_separation(np.array([1.0, 10.0]),
                                             np.array([1.0, 1.0]),
                                             np.array([1.0, 1.0]))

    assert len(a) == 2
    assert a[1] > a[0]  # longer period -> wider orbit
