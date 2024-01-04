from copy import deepcopy
import distutils
from distutils.dir_util import copy_tree
from decimal import Decimal
import multiprocessing as mp
from os import path
import pkgutil
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from tqdm import tqdm

import kaitiaki

########################################################################
#                                NOTICE                                #
#                              DEPRECATED                              #
#                              DO NOT USE                              #
########################################################################

def _allocate_cores(reserve_core: bool):
    offset = 0

    if reserve_core: offset = int(mp.cpu_count() / 4)

    return mp.cpu_count() - offset

def _worker_evolve(directory, timeout, mass):
    # our paralellised worker function. Shouldn't be necessary
    # for other uses.
    STARS = kaitiaki.STARS.STARSController(run_bs="../../..")

    out, err, reason = STARS.run(timeout=timeout,
                                 cwd=directory)

    if reason == 'finished':
        if 0.8 <= mass <= 1.6:
            STARS.evolve_through_helium_flash(timeout=timeout,
                                              filedir=directory)

    return out, err, reason

def _worker(directory, timeout):
    # our paralellised worker function. Shouldn't be necessary
    # for other uses.
    STARS = kaitiaki.STARS.STARSController(run_bs="../../..")

    out, err, reason = STARS.run(timeout=timeout,
                                 cwd=directory)

    return out, err, reason

def _worker_metal_evo(directory, timeout):
    # our paralellised worker function. Shouldn't be necessary
    # for other uses.
    STARS = kaitiaki.STARS.STARSController(run_bs="../../..")

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
        TypeError -- if ZAMS_masses is not an iterable / modin or COTable
                     location arguments are of invalid types.
        ValueError -- if ZAMS_masses is empty or metallicity is malformed.
        FileNotFoundError -- if we do not have a COtable or modin file
                             for the requested metallicity.
    """
    def __init__(self, ZAMS_masses, metallicity='z020',
                                    modin_location=None,
                                    COTable_location=None,
                                    are_logged=False):

        raise Exception("GridMaker has been deprecated. More tools will be provided in future.")
        # Most of this function is determining the validity of input.
        try:
            a = ZAMS_masses[0]
        except IndexError:
            raise TypeError("ZAMS_masses must be iterable (e.g., list, tuple).")

        if len(ZAMS_masses) == 0:
            raise ValueError("ZAMS_masses must have length >= 1.")

        # Store the masses.
        self._masses = ZAMS_masses

        Z = str(metallicity)

        self.masses_logged = are_logged

        if Z[0] == 'z':
            # We have something "like" z020. Do I have a modin for it?

            # Case I: User has not supplied a modin.
            if modin_location == None:
                try:
                    # Attempt to read a kaitiaki internal file
                    backupmodin = f"../backup_data/modins/modin.bak.{Z}"
                    pkgutil.get_data(__name__, backupmodin)
                    self._modin_location = (backupmodin, "int")
                except FileNotFoundError:
                    exc = "I don't have a modin for that metallicity."
                    raise ValueError(exc)
            # Case II: User has supplied a modin.
            elif type(modin_location) == str:
                if path.exists(modin_location):
                    # should check if properly formatted modin file
                    self._modin_location = (modin_location, "ext")
                    pass
                else:
                    exc = ("The modin file you provided ({modin_location})"
                           " does not exist.")
                    raise FileNotFoundError(exc)
            else:
                raise TypeError("The modin_location argument is invalid.")

            # do I have an opacity table for it?
            # Case I: User has not supplied a COTable.
            if COTable_location == None:
                try:
                    # Attempt to read in a kaitiaki internal file
                    backupCOtable = f"../backup_data/COtables/COtables_{Z}"
                    pkgutil.get_data(__name__, backupCOtable)
                    self._COtable_location = (backupCOtable, "int")
                except FileNotFoundError:
                    exc = "I don't have a COtable for that metallicity."
                    raise ValueError(exc)

            # Case II: User has supplied a COTable.
            elif type(COTable_location) == str:
                if path.exists(COTable_location):
                    # should check if properly formatted COTable file
                    self._COtable_location = (COTable_location, "ext")
                    pass
                else:
                    exc = ("The COTable file you provided ({COTable_location})"
                           " does not exist.")
                    raise FileNotFoundError(exc)
            else:
                raise TypeError("The COTable_location argument is invalid.")

            # Store the metallicity and the formatted one.
            self.metallicity = Z
            self.Z_fmt = float("0." + Z[1:])
        else:
            raise ValueError(("Incorrectly formatted metallicity. "
                              "I expect a format like \"z020\" for solar "
                              "(0.020)."))

        self._grid_status = 'initialized'

    def get_status(self):
        return self._grid_status

    def __str__(self):
        return f"""GridMaker Object.

Masses: {self._masses}"""

    def declare_dirname(self, dirname):
        self._models_dir = f"models/{self.metallicity}/{dirname}"

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
            {bool} -- True if make_grid() can be run, False otherwise.

        Notes:
            This creates a large amount of files. In particular,
            it creates a directory for each star to be evolved. This
            function also creates a modin and datafile in there, and
            further functions in this class such as make_grid()
            creates the files that STARS outputs.

        File Operations:
            This function sets the following parameters in data.
                NM2:     199
                IML1:    9
                RML:     <the mass of the model>
                ITH:     0
                IX:      0
                IY:      0
                IZ:      0
                ISTART:  0
                ZS:      <the grid metallicity>

            If any of these are wrong, you must edit the datafile
            yourself afterwards.
        """
        s = len(self._masses)
        kaitiaki.debug("warning", (f"I will be creating {s} subdirectories, "
                                    "located in models/."), fatal=True)

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
            STARS.terminal_command(f"mkdir models/{self.metallicity}")

            # Shuffle COTables around. We need the COtable kaitiaki stores.
            if self._COtable_location[1] == "int":
                COTable = f"../backup_data/COtables/COtables_{self.metallicity}"
                cotables = pkgutil.get_data(__name__, COTable)
                cotables = cotables.decode("utf-8")
            else:
                with open(self._COtable_location[0], 'r') as f:
                    cotables = f.read()

            # write it to file.
            with open(f"dat/COtables", "w") as f:
                f.write(cotables)

            self.declare_dirname(dirname)

            for mass in self._masses:
                # allow a wildcard * to replace the mass of the system.
                # So m* becomes m8 for mass=8.
                directory = dirname.replace('*', str(mass))

                if self.masses_logged: mass = 10**mass

                # Move modin to this subdirectory
                metal = self.metallicity

                STARS.terminal_command(f"mkdir models/{metal}/{directory}")

                if self._modin_location[1] == 'int':
                    modin_bak = f"../backup_data/modins/modin.bak.{metal}"
                    modin = pkgutil.get_data(__name__, modin_bak)
                    modin = modin.decode("utf-8")
                else:
                    with open(self._modin_location[0], 'r') as f:
                        modin = f.read()

                # Write to file.
                with open(f"models/{metal}/{directory}/modin", "w") as f:
                    f.write(modin)

                # Move data to this subdirectory. User can overwrite
                # this if they want.
                data = pkgutil.get_data(__name__, f"../backup_data/data.bak")
                data = data.decode("utf-8")

                # Write to file.
                with open(f"models/{metal}/{directory}/data", "w") as f:
                    f.write(data)

                f = f'models/{metal}/{directory}/data'

                # Set default parameters.
                # I thought about making this customisable, but
                # it's such a fundamental thing that we shan't.
                # Anyway, suppose we set NM2 to 199 when it should be
                # 499. The user can override it by editing the
                # datafile themselves - e.g. dfile.set("NM2", 499)
                with kaitiaki.file.data(f) as dfile:
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

            self._grid_status = 'setup'

            return True
        else:
            kaitiaki.debug("info", "Grid setup aborted.")

            return False

    def make_grid(self, timeout=1200, reserve_core=True):
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
            reserve_core {bool} -- Whether to reserve a core for the
                                   user to interact with their machine
                                   (True) or not (False). Set to False
                                   if you plan on running this on a
                                   cluster: it's only useful if you don't
                                   want to brick your work machine doing
                                   this code. (default: True)

        Returns:
            list -- a list corresponding to the outputs. Each element is
                    a 2-tuple, containing (mass, response), where
                    response is the output of a call to
                    STARSController.terminal_command, i.e., a 3-tuple: stdout, stderr, reason
        """

        cpu = _allocate_cores(reserve_core)
        pool = mp.Pool(cpu)

        size = len(self._masses)
        kaitiaki.debug('info', (f"Generating a grid of size {size}. "
                                f"I am using {cpu} cores. This could "
                                 "take some time!"))

        outputs = []

        outs = []

        pbar = tqdm(total=size)

        for mass in self._masses:
            # replace * with masses. Same structure as previously.
            try:
                directory = self._models_dir.replace('*', str(mass))
            except AttributeError:
                kaitiaki.debug("error", ("Model grid not initialized. "
                                         "Attempting recovery with default "
                                         "parameters..."))

                self.setup_grid()

            if self.masses_logged: mass = 10**mass

            update = lambda *a: pbar.update()
            out = pool.apply_async(_worker,
                                   args=(directory,timeout),
                                   callback=update)

            outputs.append((mass, out))

        for r in outputs:
            outs.append((r[0], r[1].get()))

        self._grid_status = 'made'

        return outs

    def evolve_grid(self, timeout=45*60, reserve_core=True):
        """Executes STARS to evolve a ZAMS grid until the end of evolution.

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
            reserve_core {bool} -- Whether to reserve a core for the
                                   user to interact with their machine
                                   (True) or not (False). Set to False
                                   if you plan on running this on a
                                   cluster: it's only useful if you don't
                                   want to brick your work machine doing
                                   this code. (default: True)

        Returns:
            list -- a list corresponding to the outputs. Each element is
                    a 2-tuple, containing (mass, response), where
                    response is the output of a call to
                    STARSController.terminal_command, i.e., a 3-tuple: stdout, stderr, reason
        """

        cpu = _allocate_cores(reserve_core)
        pool = mp.Pool(cpu)

        size = len(self._masses)
        kaitiaki.debug('info', (f"Properly Evolving a grid of size {size}. "
                                f"I am using {cpu} cores. This could "
                                 "take some time!"))

        outputs = []

        outs = []

        if self._grid_status == 'initialized':
            kaitiaki.debug("info", "I have to set up the grid first.")
            self.setup_grid()

        if self._grid_status == 'setup':
            kaitiaki.debug('info', "Grid not made! Attmepting recovery...")
            self.make_grid()

        pbar = tqdm(total=size)

        for mass in self._masses:
            # replace * with masses. Same structure as previously.
            try:
                directory = self._models_dir.replace('*', str(mass))
            except AttributeError:
                if self._grid_status == 'initialized':
                    kaitiaki.debug("error", ("Model grid not set up. "
                                             "Attempting recovery with default "
                                             "parameters..."))

                    self.setup_grid()

            if self.masses_logged: mass = 10**mass

            update = lambda *a: pbar.update()
            out = pool.apply_async(_worker_evolve,
                                   args=(directory,timeout,mass),
                                   callback=update)

            outputs.append((mass, out))

        for r in outputs:
            outs.append((r[0], r[1].get()))

        self._grid_status = 'evolved'

        return outs

    def validate_generated_grid(self):
        # Todo: implement.
        pass

    def evolve_grid_to_metallicity(self, target_metallicity,
                                         COTable_location=None,
                                         timeout=1200,
                                         reserve_core=True):
        """Evolves this grid to a target metallicity.

        TODO: document

        """

        if target_metallicity[0] != 'z':
            raise ValueError("Incorrectly formatted metallicity.")

        ZS = "0." + target_metallicity[1:]
        CH = 0.75 - 2.5 * float(ZS)

        if COTable_location == None:
            try:
                # Attempt to read in a kaitiaki internal file
                backupCOtable = f"../backup_data/COtables/COtables_{target_metallicity}"
                pkgutil.get_data(__name__, backupCOtable)
                COLocation = (backupCOtable, "int")
            except FileNotFoundError:
                exc = "I don't have a COtable for that metallicity."
                raise ValueError(exc)

        # Case II: User has supplied a COTable.
        elif type(COTable_location) == str:
            if path.exists(COTable_location):
                # should check if properly formatted COTable file
                COLocation = (COTable_location, "ext")
                pass
            else:
                exc = ("The COTable file you provided ({COTable_location})"
                       " does not exist.")
                raise FileNotFoundError(exc)
        else:
            raise TypeError("The COTable_location argument is invalid.")

        # Shuffle COTables around. We need the COtable kaitiaki stores.
        if COLocation[1] == "int":
            COTable = f"../backup_data/COtables/COtables_{target_metallicity}"
            cotables = pkgutil.get_data(__name__, COTable)
            cotables = cotables.decode("utf-8")
        else:
            with open(COLocation[0], 'r') as f:
                cotables = f.read()

        # write it to file.
        # with open(f"dat/COtables", "w") as f:
        #     f.write(cotables)

        size = len(self._masses)

        outputs = []
        outs = []
        pbar = tqdm(total=size)

        cpu = _allocate_cores(reserve_core)
        pool = mp.Pool(cpu)

        for mass in self._masses:

            try:
                directory = self._models_dir.replace('*', str(round(mass, 2)))

                new_dir = directory.replace(self.metallicity,
                                              target_metallicity)

                if not path.exists(new_dir):
                    copy_tree(directory, new_dir)

                directory = new_dir

            except AttributeError:
                raise Exception("dirname not declared.")
            except distutils.errors.DistutilsFileError as err:
                kaitiaki.debug("error", str(err))
                continue

            if self.masses_logged: mass = round(10**mass, 2)

            # Change NCH to 3
            # Change Z to 3.00E-2
            # Change CH to 0.75 - 2.5*Z
            # Load in grid
            # run_bs with COtable for z030

            datadir = f"{directory}/data"

            with open(f"{directory}/COtables", "w") as f:
                f.write(cotables)

            if not path.exists(datadir):
                data = pkgutil.get_data(__name__, f"../backup_data/data.bak")
                data = data.decode("utf-8")

                with open(f"{directory}/data", "w") as f:
                    f.write(data)

            with kaitiaki.file.data(datadir) as dfile:
                dfile.set('NCH', 3)
                dfile.set('ZS', ZS)
                dfile.set('CH', CH)
                dfile.set('ISTART', '1')
                dfile.set('IX', '1')
                dfile.set('IY', '1')
                dfile.set('IZ', '1')
                dfile.set('ITH', '0')
                dfile.set('IML1', '9')
                dfile.set('RML', mass)
                dfile.set('ISTART', '1')

            update = lambda *a: pbar.update()
            out = pool.apply_async(_worker_metal_evo,
                                   args=(directory,timeout),
                                   callback=update)

            outputs.append((mass, out))

        for r in outputs:
            outs.append((r[0], r[1].get()))

        return outs
