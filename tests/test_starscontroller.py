"""Tests for kaitiaki.STARSController.STARSController.

IMPORTANT -- read before adding to this file:
    STARSController.run() shells out to a compiled STARS executable
    (`run_bs`), which can take anywhere from seconds to hours and writes
    large files to disk. No test here (or anywhere in this suite -- see
    conftest.py's `_block_real_subprocess` fixture) may invoke it for
    real. Every test that reaches run() or terminal_command() monkeypatches
    `kaitiaki.terminal.execute` -- the single choke point both go through
    -- with a fake that returns canned output instead of touching the
    filesystem's `run_bs` binary at all. If you add a test that exercises
    a new code path through run(), mock kaitiaki.terminal.execute the same
    way; do not rely only on the global subprocess.Popen guard as your
    only safety net.
"""
import os
import shutil

import pytest

import kaitiaki
from kaitiaki.STARSController import STARSController


@pytest.fixture
def data_path(data_bak_path, tmp_path):
    """A fresh, disposable copy of the real data file at tmp_path/'data'."""
    dest = tmp_path / "data"
    shutil.copy(data_bak_path, dest)
    return dest


def test_configure_parameters_lowercases_keys():
    controller = STARSController()

    controller.configure_parameters({"NCH": 3, "iml1": 5})

    assert controller._params == {"nch": 3, "iml1": 5}


def test_configure_parameters_merges_across_calls():
    controller = STARSController()

    controller.configure_parameters({"NCH": 3})
    controller.configure_parameters({"IML1": 5})

    assert controller._params == {"nch": 3, "iml1": 5}


def test_update_datafile_and_update_run_bs():
    controller = STARSController()

    controller.update_datafile("/some/other/data")
    controller.update_run_bs("/some/other/run_bs_dir")

    assert controller._datafile == "/some/other/data"
    assert controller._run_bs_location == "/some/other/run_bs_dir"


def test_commit_parameters_writes_configured_values_and_backs_up(data_path):
    controller = STARSController()
    controller.update_datafile(str(data_path))
    controller.configure_parameters({"NCH": 8})

    controller.commit_parameters()

    backup_path = data_path.parent / (data_path.name + ".bak")
    assert backup_path.is_file()

    with kaitiaki.file.data(str(data_path)) as dfile:
        assert dfile.get("nch") == 8


def test_setup_single_evolution_writes_id_block_and_imode(data_path,
                                                           tmp_path,
                                                           monkeypatch):
    # setup_single_evolution() also calls set_period(), which works
    # against a `modin` file in the *current directory* rather than
    # wherever `dfile` lives -- so we need one on disk in cwd too.
    monkeypatch.chdir(tmp_path)
    (tmp_path / "modin").write_text("A" * 46 + "B" * 14 + "C" * 90 + "\n")

    controller = STARSController()

    controller.setup_single_evolution(dfile=str(data_path))

    lines = data_path.read_text().splitlines()
    assert lines[4].strip().startswith("6  7  0  3  0 76")
    assert controller._params["imode"] == 1


def test_setup_binary_evolution_writes_id_block_and_imode(data_path):
    # (setup_binary_evolution(), unlike setup_single_evolution(), does not
    # also call set_period(), so no modin file is needed here.)
    controller = STARSController()

    controller.setup_binary_evolution(dfile=str(data_path))

    lines = data_path.read_text().splitlines()
    assert lines[4].strip().startswith("14 14  0  9  1102")
    assert controller._params["imode"] == 2


def test_setup_evolution_rejects_unknown_mode(data_path):
    controller = STARSController()
    controller.update_datafile(str(data_path))

    with pytest.raises(ValueError):
        controller.setup_evolution(mode="not_a_real_mode")


def test_get_binary_mode(data_path, tmp_path):
    directory = tmp_path / "run"
    directory.mkdir()
    shutil.copy(data_path, directory / "data")

    controller = STARSController()

    with kaitiaki.file.data(str(directory / "data")) as dfile:
        dfile.set("IMODE", 1)
    assert controller.get_binary_mode(str(directory)) is False

    with kaitiaki.file.data(str(directory / "data")) as dfile:
        dfile.set("IMODE", 2)
    assert controller.get_binary_mode(str(directory)) is True


def test_set_period_rewrites_only_the_period_field(tmp_path):
    before = "A" * 46
    old_middle = "B" * 14
    after = "C" * 90
    modin = tmp_path / "modin"
    modin.write_text(before + old_middle + after + "\n")

    controller = STARSController()
    controller.set_period(3.5, directory=str(tmp_path),
                          boost_max_nmodels=False)

    new_line = modin.read_text().splitlines()[0]
    assert new_line[:46] == before
    assert new_line[46:48] == "  "
    assert new_line[48:60] == "{:12.6E}".format(3.5)
    assert new_line[60:] == after


def test_set_period_can_also_boost_max_nmodels(tmp_path):
    line = "A" * 200
    modin = tmp_path / "modin"
    modin.write_text(line + "\n")

    controller = STARSController()
    controller.set_period(1.0, directory=str(tmp_path),
                          boost_max_nmodels=True)

    new_line = modin.read_text().splitlines()[0]
    assert new_line[:60].endswith("{:12.6E}".format(1.0))
    assert new_line[94:94 + len(" 99999      0")] == " 99999      0"
    assert new_line[107:] == line[107:]


def test_set_period_also_updates_modin2_in_binary_mode(tmp_path):
    for name in ("modin", "modin2"):
        (tmp_path / name).write_text("A" * 46 + "B" * 14 + "C" * 90 + "\n")

    controller = STARSController()
    controller.configure_parameters({"IMODE": 2})
    controller.set_period(7.0, directory=str(tmp_path))

    for name in ("modin", "modin2"):
        new_line = (tmp_path / name).read_text().splitlines()[0]
        assert new_line[48:60] == "{:12.6E}".format(7.0)


def test_blit_writes_data_and_cotables_and_configures_zs(cotables_dir,
                                                          tmp_path):
    controller = STARSController()

    controller.blit(ZS="z001", directory=str(tmp_path))

    assert (tmp_path / "data").is_file()
    assert (tmp_path / "COtables").is_file()
    assert controller._params["zs"] == pytest.approx(0.001)
    assert controller._params["ch"] == pytest.approx(0.75 - 2.5 * 0.001)


def test_fetch_datafile_returns_shipped_data_bak(backup_data_dir):
    controller = STARSController()

    contents = controller.fetch_datafile()

    assert contents == (backup_data_dir / "data.bak").read_text()


def test_run_calls_terminal_execute_and_restores_cwd(monkeypatch,
                                                      data_path, tmp_path):
    """The one test that comes closest to run(): terminal.execute() is
    fully mocked, so run_bs is never actually invoked."""
    calls = []

    def fake_execute(command, timeout=5 * 60, cwd=None, warn=True):
        calls.append({
            "command": command,
            "timeout": timeout,
            "cwd": cwd,
            "warn": warn,
            "process_cwd": os.getcwd(),
        })
        return "out", "err", "finished"

    monkeypatch.setattr(kaitiaki.terminal, "execute", fake_execute)

    controller = STARSController(run_bs=".")
    controller.update_datafile(str(data_path))
    original_cwd = os.getcwd()

    run_dir = tmp_path / "rundir"
    run_dir.mkdir()

    stdout, stderr, reason = controller.run(cwd=str(run_dir), timeout=42)

    assert (stdout, stderr, reason) == ("out", "err", "finished")
    assert len(calls) == 1
    assert calls[0]["command"] == "./run_bs"
    assert calls[0]["timeout"] == 42
    # run() chdirs into `cwd` itself before calling terminal_command(),
    # rather than forwarding cwd through to it:
    assert calls[0]["cwd"] is None
    assert calls[0]["process_cwd"] == str(run_dir)
    # ... and always restores the original directory afterwards, even
    # though we never got as far as a real run_bs process:
    assert os.getcwd() == original_cwd


def test_run_commits_configured_parameters_before_running(monkeypatch,
                                                           data_path):
    monkeypatch.setattr(kaitiaki.terminal, "execute",
                        lambda *a, **kw: ("out", "err", "finished"))

    controller = STARSController()
    controller.update_datafile(str(data_path))
    controller.configure_parameters({"NCH": 4})

    controller.run()

    with kaitiaki.file.data(str(data_path)) as dfile:
        assert dfile.get("nch") == 4


def test_run_returns_timing_when_requested(monkeypatch, data_path):
    monkeypatch.setattr(kaitiaki.terminal, "execute",
                        lambda *a, **kw: ("out", "err", "finished"))

    controller = STARSController()
    controller.update_datafile(str(data_path))

    result = controller.run(time_me=True)

    assert len(result) == 4
    out, err, reason, delta_time = result
    assert reason == "finished"
    assert delta_time >= 0
