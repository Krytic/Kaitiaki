"""Tests for kaitiaki.kipp.cli.main.

Notes:
    kaitiaki/kipp/__main__.py runs `raise SystemExit(main())` at import
    time with no arguments, even though main() requires a `plotfile`
    argument -- so `python -m kipp` cannot actually work in the module's
    current state (the kipp integration is still being folded into
    Kaitiaki proper, per kaitiaki/kipp/NOTICE.md). We test cli.main()
    directly instead of importing __main__, since importing __main__
    would raise on import as a side effect.
"""
import numpy as np

from kaitiaki.kipp.cli import main


def _toy_data(n=4):
    return {
        "timestep": np.arange(n, dtype=float),
        "age": np.linspace(0.0, 10.0, n),
        "M": np.linspace(5.0, 4.9, n),
        "He_core": np.zeros(n),
        "CO_core": np.zeros(n),
        "conv": np.zeros((n, 12)),
        "conv_env": np.full(n, np.nan),
    }


def test_main_missing_file_returns_1_and_reports_on_stderr(tmp_path,
                                                            capsys):
    missing = tmp_path / "does_not_exist"

    exit_code = main(str(missing))

    assert exit_code == 1
    err = capsys.readouterr().err
    assert "error reading" in err


def test_main_malformed_file_returns_1_and_reports_on_stderr(tmp_path,
                                                              capsys):
    bad_plot = tmp_path / "plot"
    bad_plot.write_text("only a few tokens\n")

    exit_code = main(str(bad_plot))

    assert exit_code == 1
    err = capsys.readouterr().err
    assert "error parsing" in err

