from copy import deepcopy
from decimal import Decimal
from os import path
import time
import pkgutil
import multiprocessing as mp

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

import kaitiaki

def _worker(directory, timeout):
    # our paralellised worker function. Shouldn't be necessary
    # for other uses.
    STARS = kaitiaki.STARS.STARSController(run_bs="../..")

    out, err, reason = STARS.run(timeout=timeout,
                                 cwd=directory)

    return out, err, reason

class GridMaker:
    """
    Object to create a grid of ZAMS stars at a given metallicity.

    Arguments:
        ZAMS_Masses {iterable} -- the list of ZAMS masses we wish to
                                  generate.

    Keyword Arguments:
        metallicity {string} -- the metallicity to load. Must be
                                BPASS-formatted, e.g., z020 for solar,
                                zem5 for 10^-5.

    Raises:
        TypeError -- if ZAMS_masses is not an iterable.
        ValueError -- if ZAMS_masses is empty or metallicity is malformed.
        FileNotFoundError -- if we do not have a COtable or modin file
                             for the requested metallicity.
    """
    def __init__(self, ZAMS_masses, metallicity='z020'):
        # Most of this function is determining the validity of input.
        try:
            a = ZAMS_masses[0]
        except IndexError:
            raise TypeError("ZAMS_masses must be iterable (e.g., list, tuple).")

        if len(b) == 0:
            raise ValueError("ZAMS_masses must have length >= 1.")

        # Store the masses.
        self._masses = ZAMS_masses

        Z = str(metallicity)

        # TODO: allow user to provide own modins or cotables.
        if Z[0] == 'z':
            # We have something "like" z020. Do I have a modin for it?
            try:
                backupmodin = f"../backup_data/modins/modin.bak.{Z}"
                pkgutil.get_data(__name__, backupmodin)
            except FileNotFoundError:
                raise ValueError("I don't have a modin for that metallicity.")

            # do I have an opacity table for it?
            try:
                backupmodin = f"../backup_data/COtables/COtables_{Z}"
                pkgutil.get_data(__name__, backupmodin)
            except FileNotFoundError:
                raise ValueError("I don't have a COtable for that metallicity.")

            # Store the metallicity and the formatted one.
            self.metallicity = Z
            self.Z_fmt = float("0." + Z[1:])
        else:
            raise ValueError(("Incorrectly formatted metallicity. "
                              "I expect a format like \"z020\" for solar "
                              "(0.020)."))

    def setup_grid(self, assume_yes=False, dirname='m*'):
        """Sets up the model grid.

        Performs necessary setup for the model grid. This will create a
        lot of different model subdirectories, located in the models/
        folder. This will prompt for confirmation (requiring a "y" or
        "n" character to be input) before making filesystem changes.

        Keyword Arguments:
            assume_yes {bool} -- Whether to assume that you consent to
                                 the prompt or not. (default: False)
            dirname {string} -- The format of the directory names to
                                create. * will be replaced with the mass
                                of the star being evolved.

        Returns:
            None

        Notes:
            This creates a large amount of files. In particular,
            it creates a directory for each star to be evolved. This
            function also creates a modin and datafile in there, and
            further functions in this class such as make_grid()
            creates the files that STARS outputs.
        """
        s = len(self._masses)
        kaitiaki.debug("warning", (f"I will be creating {s} subdirectories, "
                                    "located in models/."))

        if not assume_yes: # Ask the user if we can generate all these
                           # files.
            response = None

            while response not in ['y', 'n']:
                response = input("Do you wish to continue? [y/n] ")
        else: # or ignore it if the user has indicated consent.
            kaitiaki.debug('status', "assume_yes set, bypassing consent screen")
            response = 'y'

        STARS = kaitiaki.STARS.STARSController()

        if response == 'y':
            # Create my model directory.
            STARS.terminal_command("mkdir models")

            # Shuffle COTables around. We need the COtable kaitiaki stores.
            COTable = f"../backup_data/COtables/COtables_{self.metallicity}"
            cotables = pkgutil.get_data(__name__, COTable)
            cotables = cotables.decode("utf-8")

            # write it to file.
            with open(f"dat/COtables", "w") as f:
                f.write(cotables)

            self._models_dir = f"models/{dirname}"

            for mass in self._masses:
                # allow a wildcard * to replace the mass of the system.
                # So m* becomes m8 for mass=8.
                directory = dirname.replace('*', str(mass))
                STARS.terminal_command(f"mkdir models/{directory}")

                # Move modin to this subdirectory
                metal = self.metallicity

                modin_bak = f"../backup_data/modins/modin.bak.{metal}"
                modin = pkgutil.get_data(__name__, modin_bak)
                modin = modin.decode("utf-8")

                # Write to file.
                with open(f"models/{directory}/modin", "w") as f:
                    f.write(modin)

                # Move data to this subdirectory
                data = pkgutil.get_data(__name__, f"../backup_data/data.bak")
                data = data.decode("utf-8")

                # Write to file.
                with open(f"models/{directory}/data", "w") as f:
                    f.write(data)

                f = f'models/{directory}/data'

                # Set default parameters.
                # TODO: Should this be customisable?
                with kaitiaki.datafile.DataFileParser(f) as dfile:
                    dfile.set('NM2', 199)
                    dfile.set('IML1', 9)
                    dfile.set('RML', float(mass))
                    dfile.set('ITH', 0)
                    dfile.set('IX', 0)
                    dfile.set('IY', 0)
                    dfile.set('IZ', 0)
                    dfile.set('ISTART', 0)
                    dfile.set('ZS', self.Z_fmt)

            # Now we need to edit run_bs. This in principle only needs
            # to be done once.
            pwd = STARS.terminal_command("pwd")[0]

            # read run_bs...
            with open("run_bs", "r+") as file:
                contents = []
                for l in file.readlines():
                    if l.startswith("set BSDIR"):
                        l = f"set BSDIR={pwd}\n" # redefine the BSDIR
                                                 # environment variable.
                    contents.append(l)

            # ...and write to it.
            with open("run_bs", "w") as file:
                file.seek(0)
                file.writelines(contents)

            kaitiaki.debug("info", f"Grid set up.")
        else:
            kaitiaki.debug("info", "Grid setup aborted.")

    def make_grid(self, timeout=1200):
        """Executes STARS to generate a ZAMS grid.

        Performs the actual evolution. This function is quasi-blocking:
        processes will run independently, however it shan't return until
        every instance terminates.

        This function is parallelized and can run a task on each CPU
        core available.

        If the grid has not been initialized with setup_grid(), this
        will attempt setup. If you've already set up a grid in a previous
        script, this will change nothing (you will have to indicate
        consent, this function does not assume that).

        Keyword Arguments:
            timeout {float} -- the timeout (in seconds) to apply to
                               run_bs. (default: 1200, 20 minutes.)

        Returns:
            list -- a list corresponding to the outputs. Each element is
                    a 2-tuple, containing (mass, response), where
                    response is the output of a call to
                    STARSController.terminal_command, i.e., a 3-tuple:
                        stdout, stderr, reason
        """
        cpu = mp.cpu_count()
        pool = mp.Pool(cpu)

        size = len(self._masses)
        kaitiaki.debug('info', (f"Generating a grid of size {size}. "
                                 "I am using {cpu} cores. This could "
                                 "take some time!"))

        outputs = []

        outs = []

        for mass in self._masses:
            # replace * with masses. Same structure as previously.
            try:
                directory = self._models_dir.replace('*', str(mass))
            except AttributeError:
                kaitiaki.debug("error", ("Model grid not initialized. "
                                         "Attempting recovery with default "
                                         "parameters..."))

                self.setup_grid()

            out = pool.apply_async(_worker, args=(directory,timeout))
            outputs.append((mass, out))

        for r in outputs:
            outs.append((r[0], r[1].get()))

        return outs