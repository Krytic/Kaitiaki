"""Tests for kaitiaki.kicks: supernova natal-kick velocity distributions.

The scipy-backed continuous distributions (Hobbs, MultimodalHobbs,
Verbunt) are tested for shape/sanity -- correct sample count, support
within the declared [0, 2000] km/s bounds -- rather than pinned to exact
values, since they're draws from a continuous random distribution. Sample
sizes are kept small because scipy has to numerically integrate the custom
`_pdf` to sample from these (no closed-form `_rvs`/`_ppf` is defined), and
we only need enough draws to sanity-check the plumbing here, not to
characterise the distribution.
"""
import numpy as np
import pytest

import kaitiaki


def test_sample_hobbs_returns_requested_count_and_in_support():
    np.random.seed(0)

    sample = kaitiaki.kicks.sample("Hobbs", 30)

    assert len(sample) == 30
    assert (sample >= 0).all()
    assert (sample <= 2000).all()


def test_sample_verbunt_returns_requested_count_and_in_support():
    np.random.seed(0)

    sample = kaitiaki.kicks.sample("Verbunt", 20)

    assert len(sample) == 20
    assert (sample >= 0).all()
    assert (sample <= 2000).all()


def test_sample_unknown_distribution_raises():
    # sample() first tries to build a scipy distribution class for the
    # name (raising ValueError if none exists), then falls back to
    # looking for a plain function of that name -- which raises
    # AttributeError instead if that also doesn't exist. Both are
    # "this distribution name isn't recognised" outcomes; we assert on
    # the one actually produced for a name matching neither.
    with pytest.raises(AttributeError):
        kaitiaki.kicks.sample("NotARealDistribution", 10)


def test_sample_bray_requires_mej_and_mrem():
    with pytest.raises(ValueError):
        kaitiaki.kicks.sample("Bray2018", 5)


def test_sample_bray2018_dispatches_with_given_kwargs():
    result = kaitiaki.kicks.sample("Bray2018", 3, mej=2.0, mrem=1.4)

    # mrem <= 2.5, so mns == mrem in the Bray2018 formula.
    expected = 100 * 2.0 / 1.4 - 170 * (1.4 / 1.4)
    assert len(result) == 3
    assert result == pytest.approx([expected] * 3)


def test_sample_braycustom_requires_alpha_and_beta():
    with pytest.raises(ValueError):
        kaitiaki.kicks.sample("BrayCustom", 3, mej=2.0, mrem=1.4)


def test_sample_braycustom_dispatches_with_given_kwargs():
    result = kaitiaki.kicks.sample("BrayCustom", 2, mej=2.0, mrem=1.4,
                                   alpha=10.0, beta=5.0)

    expected = 10.0 * 2.0 / 1.4 + 5.0 * (1.4 / 1.4)
    assert result == pytest.approx([expected] * 2)


def test_pdf_returns_callable_and_nonnegative():
    pdf = kaitiaki.kicks.pdf("Hobbs")

    assert callable(pdf)
    assert pdf(0.0) >= 0
    assert pdf(300.0) >= 0
