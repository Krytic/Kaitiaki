"""Tests for kaitiaki.file_handlers.datafile.DataFileParser, using a
verbatim copy of a real STARS `data` file (tests/fixtures/data.bak) as the
starting point.

Rather than hardcoding the exact column offsets baked into
kaitiaki.constants.dfile_struct (easy to get wrong by hand), most of these
tests exercise round-trip behaviour: set a parameter, then confirm get()
reflects it, on a scratch copy of the real file so the shipped data.bak
fixture is never mutated.
"""
import shutil

import pytest

import kaitiaki


@pytest.fixture
def data_path(data_bak_path, tmp_path):
    """A fresh, disposable copy of the real data file at tmp_path/'data'."""
    dest = tmp_path / "data"
    shutil.copy(data_bak_path, dest)
    return dest


def test_get_nm2_matches_known_shipped_value(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    assert dfile.get("nm2") == 199


def test_get_is_case_insensitive(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    assert dfile.get("NM2") == dfile.get("nm2")


def test_get_unknown_parameter_raises_keyerror(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    with pytest.raises(KeyError):
        dfile.get("not_a_real_parameter")


def test_get_list_of_params_returns_dict(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    result = dfile.get(["nm2", "imode"])

    assert set(result.keys()) == {"nm2", "imode"}
    assert result["nm2"] == 199


def test_set_then_get_round_trips_integer_param(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set("NCH", 7)

        assert dfile.get("nch") == 7


def test_set_then_get_round_trips_float_param(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set("RML", 12.5)

        assert dfile.get("rml") == pytest.approx(12.5)


def test_set_accepts_dict_of_params(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set({"NCH": 3, "IMODE": 2})

        assert dfile.get("nch") == 3
        assert dfile.get("imode") == 2


def test_set_with_wrong_arity_raises_typeerror(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        with pytest.raises(TypeError):
            dfile.set("a", "b", "c")
        with pytest.raises(TypeError):
            dfile.set(5)


def test_set_zams_mass_sets_rml_and_iml1(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set_zams_mass(15.0)

        assert dfile.get("rml") == pytest.approx(15.0)
        assert dfile.get("iml1") == 9


def test_set_zs_from_zprefixed_string(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set("ZS", "z020")

        assert dfile.get("zs") == pytest.approx(0.020, abs=1e-4)


def test_write_persists_across_reopen(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set("NCH", 9)

    dfile2 = kaitiaki.file.data(str(data_path))

    assert dfile2.get("nch") == 9


def test_as_dict_covers_every_known_parameter(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    result = dfile.as_dict()

    assert set(result.keys()) == set(kaitiaki.constants.dfile_struct.keys())


def test_equality_operator(data_path, tmp_path):
    dfile_a = kaitiaki.file.data(str(data_path))

    other_path = tmp_path / "data_copy"
    shutil.copy(data_path, other_path)
    dfile_b = kaitiaki.file.data(str(other_path))

    assert dfile_a == dfile_b

    with kaitiaki.file.data(str(other_path)) as dfile:
        dfile.set("NCH", dfile.get("nch") + 1)

    assert dfile_a != kaitiaki.file.data(str(other_path))


def test_compare_reports_only_mismatches(data_path, tmp_path):
    other_path = tmp_path / "data_copy"
    shutil.copy(data_path, other_path)

    with kaitiaki.file.data(str(other_path)) as dfile:
        dfile.set("NCH", 42)

    dfile_a = kaitiaki.file.data(str(data_path))
    dfile_b = kaitiaki.file.data(str(other_path))

    mismatches = dfile_a.compare(dfile_b)

    assert list(mismatches.keys()) == ["nch"]


def test_backup_if_not_exists_creates_backup_once(data_path):
    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.backup_if_not_exists()

    backup_path = data_path.parent / (data_path.name + ".bak")
    assert backup_path.is_file()
    first_backup_contents = backup_path.read_text()

    with kaitiaki.file.data(str(data_path)) as dfile:
        dfile.set("NCH", 123)
        dfile.backup_if_not_exists()  # should be a no-op: backup exists

    assert backup_path.read_text() == first_backup_contents


def test_str_returns_exact_file_contents(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    assert str(dfile) == data_path.read_text()


def test_explain_known_disambiguable_parameter(data_path, capsys):
    dfile = kaitiaki.file.data(str(data_path))

    dfile.explain("cepr")

    out = capsys.readouterr().out
    assert "Common Envelope Prescription" in out


def test_show_in_file_highlights_requested_parameter(data_path, capsys):
    dfile = kaitiaki.file.data(str(data_path))

    dfile.show_in_file("nch")

    out = capsys.readouterr().out
    assert out != ""


def test_show_in_file_rejects_unknown_parameter(data_path):
    dfile = kaitiaki.file.data(str(data_path))

    with pytest.raises(KeyError):
        dfile.show_in_file("not_a_real_parameter")
