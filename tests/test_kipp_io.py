"""Tests for kaitiaki.kipp.io.load_plot.

kipp is a small, mostly self-contained module (numpy + matplotlib only,
per its own README) that decodes STARS `plot` files into shaded Kippenhahn
diagrams. These tests build minimal, synthetic plot-file lines rather than
shipping a real (multi-megabyte) plot file, since load_plot only cares
about whitespace-separated field *positions*, not fixed-width columns.
"""
import numpy as np
import pytest

from kaitiaki.kipp.io import load_plot


def _row(timestep, age, M, He_core, CO_core, conv, n_cols=23):
    """Builds one whitespace-separated STARS plot-file row.

    Args:
        timestep (int): model/timestep number (column 0).
        age (float): stellar age (column 1).
        M (float): total mass (column 5).
        He_core (float): helium core mass (column 6).
        CO_core (float): CO core mass (column 7).
        conv (list[float]): the 12 M_conv1..12 boundary values
            (columns 11..22).
        n_cols (int): total number of columns to pad the row out to with
            zeros (must be >= 23).

    Returns:
        str: a single plot-file line, newline-terminated.
    """
    assert len(conv) == 12
    fields = [0.0] * n_cols
    fields[0] = timestep
    fields[1] = age
    fields[5] = M
    fields[6] = He_core
    fields[7] = CO_core
    fields[11:23] = conv
    return " ".join(str(v) for v in fields) + "\n"


def test_load_plot_reads_core_columns(tmp_path):
    conv = [0.0] * 12
    conv[0] = -8.0  # one conv/semi boundary
    lines = [
        _row(1, 100.0, 20.0, 0.0, 0.0, conv),
        _row(2, 200.0, 19.9, 1.0, 0.5, conv),
    ]
    plotfile = tmp_path / "plot"
    plotfile.write_text("".join(lines))

    data = load_plot(plotfile)

    assert data["timestep"].tolist() == [1, 2]
    assert data["age"].tolist() == [100.0, 200.0]
    assert data["M"].tolist() == [20.0, 19.9]
    assert data["He_core"].tolist() == [0.0, 1.0]
    assert data["CO_core"].tolist() == [0.0, 0.5]
    assert data["conv"].shape == (2, 12)
    assert np.isnan(data["conv_env"]).all()  # rows too short to have it


def test_load_plot_skips_blank_lines(tmp_path):
    conv = [0.0] * 12
    lines = [_row(1, 1.0, 1.0, 0.0, 0.0, conv), "\n", "   \n",
             _row(2, 2.0, 2.0, 0.0, 0.0, conv)]
    plotfile = tmp_path / "plot"
    plotfile.write_text("".join(lines))

    data = load_plot(plotfile)

    assert len(data["timestep"]) == 2


def test_load_plot_raises_on_short_line(tmp_path):
    plotfile = tmp_path / "plot"
    plotfile.write_text("1 2 3\n")

    with pytest.raises(ValueError, match="line 1"):
        load_plot(plotfile)


def test_load_plot_recovers_merged_conv_env_token(tmp_path):
    conv = [0.0] * 12
    fields = [0.0] * 71
    fields[0] = 1
    fields[1] = 10.0
    fields[5] = 20.0
    fields[11:23] = conv
    # Simulate the fixed-width overflow described in the module docstring:
    # M_conv-env (5dp) glued to the next field with no separating space.
    fields[70] = "19.83829100.25554"
    line = " ".join(str(v) for v in fields) + "\n"

    plotfile = tmp_path / "plot"
    plotfile.write_text(line)

    data = load_plot(plotfile)

    assert data["conv_env"][0] == pytest.approx(19.83829)


def test_load_plot_missing_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        load_plot(tmp_path / "does_not_exist")
