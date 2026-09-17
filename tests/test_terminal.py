"""Tests for kaitiaki.terminal: the subprocess wrapper that
STARSController.run() (and everything else that shells out) ultimately
goes through.

None of these tests spawn a real process. subprocess.Popen is replaced
with a small fake that mimics just enough of the real interface (a context
manager with .communicate()/.wait()/.kill()/.pid) for
_custom_subprocess_handler to exercise its real control flow, including
the custom timeout-and-process-group-kill logic described in its own
docstring.
"""
import os
import subprocess

import pytest

import kaitiaki


class _FakePopen:
    """A stand-in for subprocess.Popen that never touches a real process.

    Args:
        effects (list): a queue of values consumed in order by successive
            calls to .communicate(). An entry that is a BaseException
            instance is raised instead of returned.
        pid (int): the fake pid reported via self.pid.
    """
    def __init__(self, effects, pid=99999):
        self._effects = list(effects)
        self.pid = pid
        self.killed = False
        self.waited = False

    def __call__(self, *args, **kwargs):
        """Makes the instance usable directly as subprocess.Popen."""
        self.call_args = args
        self.call_kwargs = kwargs
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def communicate(self, input=None, timeout=None):
        if self._effects:
            effect = self._effects.pop(0)
            if isinstance(effect, BaseException):
                raise effect
            return effect
        return (b"stdout", b"stderr")

    def wait(self):
        self.waited = True

    def kill(self):
        self.killed = True


@pytest.fixture
def fake_popen(monkeypatch):
    """Installs a _FakePopen instance in place of subprocess.Popen.

    Returns:
        callable: a factory `make_fake(effects=())` that installs and
        returns the fake, so a test can both configure its behaviour and
        later inspect .killed / .waited / .call_args on it.
    """
    def make_fake(effects=()):
        fake = _FakePopen(effects)
        monkeypatch.setattr(subprocess, "Popen", fake)
        return fake

    return make_fake


def test_execute_normal_completion(fake_popen):
    fake_popen(effects=[(b"hello\n", b"")])

    stdout, stderr, reason = kaitiaki.terminal.execute("echo hello")

    assert stdout == "hello"
    assert stderr == ""
    assert reason == "finished"


def test_execute_passes_cwd_and_restores_it(fake_popen, tmp_path):
    fake_popen(effects=[(b"", b"")])
    original_cwd = os.getcwd()

    kaitiaki.terminal.execute("echo hi", cwd=str(tmp_path))

    assert os.getcwd() == original_cwd


def test_execute_timeout_kills_process_group(fake_popen, monkeypatch):
    killed_pgids = []
    monkeypatch.setattr(kaitiaki.terminal.os, "getpgid", lambda pid: pid)
    monkeypatch.setattr(kaitiaki.terminal.os, "killpg",
                        lambda pgid, sig: killed_pgids.append((pgid, sig)))

    fake = fake_popen(effects=[
        subprocess.TimeoutExpired(cmd="run_bs", timeout=1),
        (b"partial", b"timed out"),
    ])

    stdout, stderr, reason = kaitiaki.terminal.execute("run_bs",
                                                        timeout=1,
                                                        warn=False)

    assert reason == "timeout"
    assert fake.waited is True
    assert killed_pgids  # os.killpg was invoked to reap the process group


def test_execute_other_exception_kills_and_reraises(fake_popen):
    fake = fake_popen(effects=[KeyboardInterrupt()])

    with pytest.raises(KeyboardInterrupt):
        kaitiaki.terminal.execute("run_bs", warn=False)

    assert fake.killed is True


def test_execute_warns_on_long_timeout(fake_popen, monkeypatch):
    fake_popen(effects=[(b"", b"")])
    messages = []
    monkeypatch.setattr(
        kaitiaki, "debug",
        lambda msgtype, message: messages.append((msgtype, message))
    )

    kaitiaki.terminal.execute("echo hi", timeout=21 * 60, warn=True)

    assert any(m[0] == "warning" for m in messages)


def test_execute_no_warning_below_threshold(fake_popen, monkeypatch):
    fake_popen(effects=[(b"", b"")])
    messages = []
    monkeypatch.setattr(
        kaitiaki, "debug",
        lambda msgtype, message: messages.append((msgtype, message))
    )

    kaitiaki.terminal.execute("echo hi", timeout=60, warn=True)

    assert messages == []
