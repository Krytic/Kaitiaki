from copy import deepcopy
from decimal import Decimal
from os import path

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

import kaitiaki

class STARSController:
    def __init__(self, verbose_output=True):
        """
        Represents an instance of the STARS code
        """
        self._verbose_output = verbose_output
        self._output_dir = "."

    def output(self, type, message):
        """Outputs a message to stdout, if verbose_output is on. Does
        nothing otherwise.

        Arguments:
            type {string} -- The message type (warning, error, info, status) to output
            message {string} -- The message to output
        """
        if type == 'warning':
            header = "\033[93m\033[1m[WARNING]\033[0m "
        elif type == 'error':
            header = "\033[91m\033[1m[ERROR]\033[0m "
        elif type == 'info':
            header = "\033[96m\033[1m[INFO]\033[0m "
        elif type == 'status':
            header = "\033[92m\033[1m[PROGRAM STATUS]\033[0m "

        if self._verbose_output:
            print(header + str(message))

    def configure_parameters(self, params):
        """Sets up the parameters. Doesn't commit them yet.

        Arguments:
            params {dict} -- A dictionary containing key: value pairs
            of what to write to the datafile.
        """
        self._params = params

    def set_output_directory(self, directory):
        """Sets up the output directory.

        If the directory doesn't exist, this will create it.

        Arguments:
            directory {string} -- The directory to output files to.
        """
        self._output_dir = directory

        if not path.exists(self._output_dir):
            self.terminal_command(f'mkdir {self._output_dir}')

    def terminal_command(self, command, timeout=5*60):
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

        proc = None

        try:
            proc = subprocess.run(cmd, timeout=timeout, capture_output=True)
        except subprocess.TimeoutExpired:
            self.output('error', f'Command Timeout. The command I was running was: {command}')
            reason = "timeout"
            stdout = b'error'
            stderr = b'error'

        if proc != None:
            # Finished!

            stdout = proc.stdout
            stderr = proc.stderr
            reason = "finished"

        return stdout.decode('utf-8').strip(), stderr.decode('utf-8').strip(), reason

    def run(self):
        """Runs ./run_bs

        Returns:
            tuple -- The output from run() above.
        """
        with kaitiaki.datafile.DataFileParser("data") as dfile:
            dfile.backup_if_not_exists()

            self.output('status', 'Writing to data file')

            for param, val in self._params.items():
                dfile.set(param, val)

        return self.terminal_command(f'./run_bs')

    def run_default_evolution(self, zams_mass):
        """Performs one run from pre-ZAMS until the end of evolution.

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
                           'ISTOP': 0}

        self.output('status', 'Configuring Parameters')
        self.configure_parameters(pre_zams_params)
        self.output('status', 'Beginning ZAMS evolution')
        self.run()

        out = self._output_dir
        self.output('status', 'Moving output files')
        self.terminal_command(f'mv plot {out}/plot.zams')
        self.terminal_command(f'mv out {out}/out.zams')
        self.terminal_command(f'tail -399 modout > modin')
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
        self.run()
        self.output('status', 'Moving output files')

        self.terminal_command(f'mv plot {out}/plot.model')
        self.terminal_command(f'mv out {out}/out.model')
        self.terminal_command(f'mv modout {out}/modout.model')
