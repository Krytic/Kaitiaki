"""Shared pytest fixtures and safety nets for the Kaitiaki test suite.

Notes:
    The most important fixture in this file is `_block_real_subprocess`.
    STARSController.run() (and everything that calls it) ultimately shells
    out to `run_bs`, a compiled STARS executable. We never want a test to
    accidentally invoke that executable -- it can take minutes to hours to
    run and writes large files to disk. `_block_real_subprocess` is
    autouse and replaces `subprocess.Popen` with a stub that raises
    immediately, for every test in the suite. Tests that need to exercise
    subprocess-adjacent code (test_terminal.py) install their own local
    monkeypatch of subprocess.Popen on top of this one, so real processes
    are never spawned there either -- see test_terminal.py's module
    docstring for details.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # never try to open a GUI window during tests

import pytest


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
FIXTURES_DIR = TESTS_DIR / "fixtures"
BACKUP_DATA_DIR = REPO_ROOT / "backup_data"


class RealSubprocessBlocked(RuntimeError):
    """Raised if a test attempts to spawn a real subprocess.

    Notes:
        This should never happen -- it means a test (or the code it is
        exercising) tried to shell out for real instead of going through a
        mock. See the module docstring for why this is blocked globally.
    """


@pytest.fixture(autouse=True)
def _block_real_subprocess(monkeypatch):
    """Prevents any test from spawning a real subprocess by default.

    Args:
        monkeypatch (pytest.MonkeyPatch): standard pytest fixture.

    Notes:
        kaitiaki.STARSController.run() shells out to the compiled STARS
        executable (`run_bs`). No test in this suite should ever do that
        for real. Individual tests that need to exercise subprocess
        plumbing (see test_terminal.py) install their own fake
        subprocess.Popen on top of this fixture from within the test
        itself; pytest/monkeypatch tears both down automatically at the
        end of the test, in last-in-first-out order.
    """
    import subprocess as _subprocess

    def _guard(*args, **kwargs):
        raise RealSubprocessBlocked(
            "A test tried to spawn a real subprocess via subprocess.Popen "
            f"(args={args!r}, kwargs={kwargs!r}). Mock subprocess.Popen or "
            "kaitiaki.terminal.execute instead of letting this through."
        )

    monkeypatch.setattr(_subprocess, "Popen", _guard)


@pytest.fixture(autouse=True)
def _close_all_figures():
    """Closes every matplotlib figure after each test.

    Notes:
        Several kaitiaki modules (STARSController, kipp.render, the
        file_handlers plot classes, ...) create matplotlib figures without
        closing them. Left unchecked across a whole test session this
        leaks memory and can eventually trigger matplotlib's "too many
        open figures" warning. This fixture keeps every test isolated.
    """
    yield
    import matplotlib.pyplot as plt
    plt.close("all")


@pytest.fixture
def repo_root():
    """Returns the repository root as a pathlib.Path."""
    return REPO_ROOT


@pytest.fixture
def backup_data_dir():
    """Returns the path to backup_data/, skipping the test if absent.

    Notes:
        backup_data/ ships with the repository but isn't guaranteed to be
        present in every checkout (e.g. a shallow clone, or a source
        distribution that excludes large data files). Tests relying on it
        should request this fixture rather than hardcoding the path.
    """
    if not BACKUP_DATA_DIR.is_dir():
        pytest.skip("backup_data/ is not present in this checkout")
    return BACKUP_DATA_DIR


@pytest.fixture
def cotables_dir(backup_data_dir):
    """Returns backup_data/COtables/, skipping the test if absent.

    Notes:
        COtables/ is ~240MB across its 25 metallicity files, so it's a
        plausible thing to be missing from a shallow clone or a stripped
        checkout even when the rest of backup_data/ is present. Tests
        that specifically need it (blit(), list_internal_files()) should
        request this fixture rather than the more general
        `backup_data_dir`.
    """
    path = backup_data_dir / "COtables"
    if not path.is_dir():
        pytest.skip("backup_data/COtables/ is not present in this checkout")
    return path


@pytest.fixture
def data_bak_path():
    """Path to a bundled, verbatim copy of a real STARS `data` file.

    Notes:
        This is a byte-for-byte copy of backup_data/data.bak, checked into
        tests/fixtures so datafile tests don't depend on backup_data/
        being present, and never risk mutating the shipped copy (tests
        that need to write to it should copy it again into tmp_path; see
        the `data_path` fixture in test_datafile.py and
        test_starscontroller.py).
    """
    path = FIXTURES_DIR / "data.bak"
    if not path.is_file():
        pytest.skip("tests/fixtures/data.bak is missing")
    return path
