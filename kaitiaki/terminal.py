import os
import time
import subprocess

import kaitiaki


def execute(command, timeout=5*60, cwd=None, warn=True):
    """Executes a terminal command.

    Has a default timeout period of 5 minutes.

    Arguments:
        command {string} -- The command to execute.

    Keyword Arguments:
        timeout {int} -- The time (in seconds) that the command
        should execute for at most. (default: {5*60})

    Returns:
        tuple -- A 3-tuple representing (stdout, stderr, reason for termination)
    """
    cmd = command.split(" ")

    timer = time.strftime('%Hh %Mm %Ss', time.gmtime(timeout))

    if timeout > 20 * 60:
        from datetime import datetime, timedelta
        now = datetime.now()
        delta = timedelta(seconds=timeout)

        TAT = now + delta

        now = now.strftime("%d %b, %H:%M:%S (%p)")
        TAT = TAT.strftime("%d %b, %H:%M:%S (%p)")

        if warn:
            err_msg = (f"I've been asked to run {command} with a maximum "
                       f"timeout of {timer} (default is 00h 05m 00s).\n"
                       f"I'm assuming this is right, but double check if you "
                       f"were not expecting this.\nCurrent time: {now}\n"
                       f"Timeout at: {TAT}")

            kaitiaki.debug("warning", err_msg)

    stdout, stderr, reason = _custom_subprocess_handler(command, timeout, cwd)

    return stdout, stderr, reason


def _custom_subprocess_handler(command, timeout=5*60, cwd=None):
    """
        Okay, this one deserves an explanation.

        I noticed that for jobs that took longer than timeout seconds
        to run, they wouldn't be killed correctly.

        https://stackoverflow.com/questions/70587181/terminate-child-process-on-subprocess-timeoutexpired/72135833

        So what I ended up finding out is that Python treats subprocesses
        which spawn groups weirdly. STARS, I suppose, does such
        a thing. Essentially, the SIGKILL signal gets sent to the
        subprocess -- *but not the group*. In turn, this means that
        the subprocesses that do spawn groups don't get killed when
        a TimeoutException is thrown. So, what follows is a copy
        of the Python source code, with the group-killer added
        (see "the magic line!" below).

        &copy; original author. Modified on 06 May 2022 by Sean Richards
        (@Krytic).
    """
    with subprocess.Popen(command.split(" "),
                          preexec_fn=os.setsid,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE) as process:

        wd = os.getcwd()

        try:
            if cwd is not None:
                # Man fuck linux
                for d in cwd.split("/"):
                    os.chdir(d)
            stdout, stderr = process.communicate(None, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            import signal

            # The magic line!
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)

            try:
                import msvcrt
            except ModuleNotFoundError:
                _mswindows = False
            else:
                _mswindows = True

            if _mswindows:
                # Windows accumulates the output in a single blocking
                # read() call run on child threads, with the timeout
                # being done in a join() on those threads.  communicate()
                # _after_ kill() is required to collect that and add it
                # to the exception.
                exc.stdout, exc.stderr = process.communicate()
            else:
                # POSIX _communicate already populated the output so
                # far into the TimeoutExpired exception.
                process.wait()
            reason = 'timeout'
            stdout, stderr = process.communicate()
        except:  # Including KeyboardInterrupt, communicate handled that.
            process.kill()
            # We don't call process.wait() as .__exit__ does that for us.
            reason = 'other'
            stdout, stderr = process.communicate()
            raise
        else:
            reason = 'finished'
        finally:
            os.chdir(wd)

        try:
            return (stdout.decode('utf-8').strip(),
                    stderr.decode('utf-8').strip(),
                    reason)
        except AttributeError:
            try:
                return stdout.strip(), stderr.strip(), reason
            except AttributeError:
                return stdout, stderr, reason

        return stdout, stderr, reason