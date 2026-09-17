"""Tests for the loose helper functions defined directly in
kaitiaki/__init__.py: format_metallicity, debug, load_file, and
list_internal_files.
"""
import pytest

import kaitiaki


@pytest.mark.parametrize("z,expected", [
    (0.020, "z020"),
    (0.001, "z001"),
    (0.0001, "zem4"),
    (0.00001, "zem5"),
])
def test_format_metallicity_from_float(z, expected):
    assert kaitiaki.format_metallicity(z) == expected


def test_format_metallicity_from_string_zprefixed_is_a_no_op():
    assert kaitiaki.format_metallicity("z020") == "z020"


def test_format_metallicity_from_plain_number_string():
    assert kaitiaki.format_metallicity("0.020") == "z020"


def test_format_metallicity_rejects_other_types():
    with pytest.raises(ValueError):
        kaitiaki.format_metallicity([0.02])


def test_debug_rejects_unknown_msgtype():
    with pytest.raises(AssertionError):
        kaitiaki.debug("not_a_real_type", "hello")


@pytest.mark.parametrize("msgtype,logger_method", [
    ("info", "info"),
    ("status", "status"),
    ("warning", "warn"),
    ("error", "error"),
])
def test_debug_dispatches_to_the_right_logger_method(monkeypatch, msgtype,
                                                      logger_method):
    calls = []
    monkeypatch.setattr(kaitiaki.log, logger_method,
                        lambda message: calls.append(message))

    kaitiaki.debug(msgtype, "hello world")

    assert calls == ["hello world"]


def test_list_internal_files_includes_expected_entries(cotables_dir):
    files = kaitiaki.list_internal_files()

    assert "data.bak" in files
    assert any(f.startswith("COtables" + "/") for f in files)


def test_load_file_returns_data_bak_contents(backup_data_dir):
    contents = kaitiaki.load_file("data.bak")

    assert contents == (backup_data_dir / "data.bak").read_text()


def test_load_file_rejects_path_traversal(backup_data_dir):
    with pytest.raises(ValueError):
        kaitiaki.load_file("../etc/passwd")


def test_load_file_rejects_disallowed_folder(backup_data_dir):
    with pytest.raises(ValueError):
        kaitiaki.load_file("not_a_real_folder/whatever")


def test_load_file_missing_file_raises_ioerror(backup_data_dir):
    with pytest.raises(IOError):
        kaitiaki.load_file("data.bak.does_not_exist_honest")
