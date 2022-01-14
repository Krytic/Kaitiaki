from copy import deepcopy
from decimal import Decimal
from os import path
import os
import time

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

import kaitiaki

class STARSController:
    def __init__(self, verbose_output=True, run_bs="."):
        """
        Represents an instance of the STARS code
        """
        self._verbose_output = verbose_output
        self._output_dir = "."
        self._params = dict()
        self._run_bs_location = run_bs
        self._datafile = "data"

    def shit_location(self, cwd):
        wd = os.getcwd()

        for d in cwd.split("/"):
            os.chdir(d)

        resp = self.terminal_command('pwd')[0]

        os.chdir(wd)

        return resp

    def output(self, type, message):
        """Outputs a message to stdout, if verbose_output is on. Does
        nothing otherwise.

        Arguments:
            type {string} -- The message type (warning, error, info, status) to output
            message {string} -- The message to output
        """
        if self._verbose_output:
            kaitiaki.debug(type, str(message))

    def configure_parameters(self, params):
        """Sets up the parameters. Doesn't commit them yet.

        Arguments:
            params {dict} -- A dictionary containing key: value pairs
            of what to write to the datafile.
        """
        for k, v in params.items():
            self._params[k.lower()] = v

    def _write_orbital_equations(self, block):
        with open(self._datafile, 'r+') as f:
            dfile = f.readlines()

            f.seek(0)

            blk = block.split("\n")
            for i in range(4):
                dfile[4+i] = blk[i] + "\n"

            f.writelines(dfile)
            f.truncate()

    def setup_binary_evolution(self):
        binary_block = """ 14 14  0  9  1102  0  0  0 99
  1  2  4 16 17 19 13 14 29  5  3  9 10 11 12 15 20 18 24 25 26 27 30  8  7  6 23 22 21  0
  7  8  9 10 11 12 14 22 23 24 25 26 27 29  4  2  1  3  5  6 19 17 16 13 28 18 20 21  0  0
  4  5  6  7  8  9 10 19 20 21 22 23 24 25  2  3  1 17 18 16  2  3  1  4  5 17 18 16 20  0"""

        self._write_orbital_equations(binary_block)

        params = {
            'imode': 2,
            'IML1': 5,
            'IML2': 5,
            'RML': 0,
            'ITH': 1,
            'IX': 1,
            'IY': 1,
            'IZ': 1,
            'ISTART': 1
        }

        self.configure_parameters(params)

    def setup_single_evolution(self):
        single_block = """  6  7  0  3  0 76  0  0  0 99
  1  2  4  5  3  9 10 11 12 15  8  7  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  7  8  9 10 11 12 14  4  2  1  3  5  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  4  5  6  7  8  9 10  2  3  1  2  3  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"""

        self._write_orbital_equations(single_block)

        params = {
            'imode': 1,
            'IML1': 5,
            'RML': 0,
            'ITH': 1,
            'IX': 1,
            'IY': 1,
            'IZ': 1,
            'ISTART': 1
        }

        self.configure_parameters(params)

    def set_output_directory(self, directory):
        """Sets up the output directory.

        If the directory doesn't exist, this will create it.

        Arguments:
            directory {string} -- The directory to output files to.
        """
        self._output_dir = directory

        if not path.exists(self._output_dir):
            self.terminal_command(f'mkdir {self._output_dir}')

    def terminal_command(self, command, timeout=5*60, cwd=None):
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

        if timeout > 10 * 60:
            self.output("warning", f"Woah! I've been asked to run {command} with a maximum timeout of {timer} (default is 00h 05m 00s). I'm assuming this is right, but double check if you were not expecting this.")

        proc = None

        wd = os.getcwd()
        try:
            if cwd is not None:
                # Man fuck linux
                for d in cwd.split("/"):
                    os.chdir(d)

            proc = subprocess.run(cmd, timeout=timeout, capture_output=True)
        except subprocess.TimeoutExpired:
            self.output('error', f'Command timed out (waited {timer}). The command I was running was: {command}')
            # proc.terminate() # can't kill children without this!
            reason = "timeout"
            stdout = b'error'
            stderr = b'error'
        finally:
            os.chdir(wd)

        if proc != None:
            # Finished!

            stdout = proc.stdout
            stderr = proc.stderr
            reason = "finished"

        return stdout.decode('utf-8').strip(), stderr.decode('utf-8').strip(), reason

    def run(self, timeout=5*60, cwd=None):
        """Runs ./run_bs

        Returns:
            tuple -- The output from run() above.
        """
        wd = os.getcwd()

        if cwd is not None:
            for d in cwd.split("/"):
                os.chdir(d)

        with kaitiaki.datafile.DataFileParser(self._datafile) as dfile:
            dfile.backup_if_not_exists()

            for param, val in self._params.items():
                dfile.set(param, val)

        out, err, reason = self.terminal_command(f'{self._run_bs_location}/run_bs', timeout=timeout)

        os.chdir(wd)

        return out, err, reason

    def modout_to_modin(self, modout_location="modout", modin_location="modin"):
        # order of preference:
        # - existing datafile
        # - declared new parameters
        # - the STARS default (199)
        if path.exists(self._datafile):
            with kaitiaki.datafile.DataFileParser(self._datafile) as dfile:
                nmesh = dfile.get('NM2')
        elif 'nm2' in self._params.keys():
            nmesh = self._params['nm2']
        else:
            nmesh = 199 # assume default

        lines = nmesh * 2 + 1

        with open(modout_location, 'r') as f:
            output = f.readlines()

        N = len(output)

        output = output[N-lines:]

        with open(modin_location, 'w') as f:
            f.writelines(output)

    def run_default_evolution(self, zams_mass):
        """Performs one run from pre-ZAMS until the end of evolution.
        Assumes NORMAL file structure.

        Arguments:
            zams_mass {float} -- The ZAMS mass we should target.
        """
        self.output('status', 'Preparing for ZAMS evolution')
        self.terminal_command("rm modin")
        self.terminal_command("rm data")

        self.terminal_command("cp modin.bak modin")
        self.terminal_command("cp data.bak data")

        pre_zams_params = {'IML1': 9,
                           'RML': float(zams_mass),
                           'ITH': 0,
                           'IX': 0,
                           'IY': 0,
                           'IZ': 0,
                           'ISTART': 0}

        self.output('status', 'Configuring Parameters')
        self.configure_parameters(pre_zams_params)
        self.output('status', 'Beginning ZAMS evolution')
        output, error, status = self.run()

        if self._verbose_output:
            self.output("status", output)
            self.output("error", error)

        out = self._output_dir
        self.output('status', 'Moving output files')
        self.terminal_command(f'mv plot {out}/plot.zams')
        self.terminal_command(f'mv out {out}/out.zams')
        self.modout_to_modin()
        self.terminal_command(f'mv modout {out}/modout.zams')

        main_evo_params = {'IML1': 5,
                           'RML': 0,
                           'ITH': 1,
                           'IX': 1,
                           'IY': 1,
                           'IZ': 1,
                           'ISTART': 1}

        self.output('status', 'Configuring Parameters')
        self.configure_parameters(main_evo_params)

        self.output('status', 'Beginning post-ZAMS evolution')

        output, error, status = self.run()

        if self._verbose_output:
            self.output("status", output)
            self.output("error", error)

        self.output('status', 'Moving output files')

        self.terminal_command(f'mv plot {out}/plot.model')
        self.terminal_command(f'mv out {out}/out.model')
        self.terminal_command(f'mv modout {out}/modout.model')
