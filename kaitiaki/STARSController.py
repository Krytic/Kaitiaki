from copy import deepcopy
from decimal import Decimal
import threading

import multiprocessing as mp
import os
import pkgutil
import requests
import shutil
import time
import traceback
import zipfile

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from tqdm import tqdm

import kaitiaki

import importlib.resources as pkg_resources

# Auxilary functions.


class ServerError(Exception):
    """Represents an error from the server

    Syntactic sugar for a server error.
    """
    def __init__(self, server_code, message):
        self.server_code = server_code
        self.message = message

        super().__init__(f'{self.server_code} {self.message}')


def install(install_path):
    """Installs STARS

    Downloads the latest version of AotearoaSTARS from GitHub.

    Args:
        install_path (str): The location to install AotearoaSTARS to.

    Raises:
        FileNotFoundError: If the repository could not be found
        ServerError: If the HTTP response code is not 200
    """
    repo = 'UoA-Stars-And-Supernovae/STARS:master'

    repo_name = repo.split('/')[1].split(":")[0]
    branch = repo.split('/')[1].split(":")[1]

    url = f'https://github.com/{repo}/archive/refs/heads/{branch}.zip'

    try:
        response = requests.get(url)
    except OSError as exc:
        # Certificate not found
        raise FileNotFoundError(f"Could not access the {repo} repo") from exc
        # I dislike it raising an OSError.

    if response.status_code == 200:
        open(f'{install_path}/{repo_name}.zip', 'wb').write(response.content)
    else:
        raise ServerError(response.status_code, response.reason)

    with zipfile.ZipFile(f'{install_path}/{repo_name}.zip', 'r') as zip_ref:
        zip_ref.extractall()

    source = f"{install_path}/{repo_name}-{branch}"

    files_list = os.listdir(source)

    for files in files_list:
        shutil.move(f"{source}/{files}", install_path)

    if os.path.realpath(source) != '/':
        # Look, this should never happen unless the user does
        # some wack shit, but I'm paranoid.
        shutil.rmtree(source)

    os.remove(f'{install_path}/{repo_name}.zip')

    kaitiaki.terminal.execute('make', cwd=install_path)


def recompile(install_path):
    """Recompiles the STARS code

    Recompiles the STARS code. Executes make clean && make to do so.

    Args:
        install_path (str): The location that the STARS code is installed to

    Returns:
        tuple: a 2-tuple of 3-tuples representing ((stdout, stderr, termination_reason), (stdout, stderr, termination_reason)) for (make clean, make) respectively.

    Raises:
        ChildProcessError: If either command fails.
    """
    o, e, r = kaitiaki.terminal.execute('make clean', cwd=install_path)

    if r == 'finished':
        o2, e2, r2 = kaitiaki.terminal.execute('make', cwd=install_path)

        if r2 == 'finished':
            return ((o, e, r), (o2, e2, r2))
        else:
            raise ChildProcessError('make failed.')
    else:
        raise ChildProcessError('make clean failed.')


def _allocate_cores(reserve_core: bool):
    offset = 0

    n_cores = mp.cpu_count()

    if reserve_core:
        if n_cores > 12:
            offset = 2
        else:
            offset = int(n_cores // 4)

    return n_cores - offset


def _worker_evolve(directory, timeout, mass, do_he_flash, run_bs, STARS):
    # our paralellised worker function. Shouldn't be necessary
    # for other uses.

    STARS.update_run_bs(run_bs)

    if os.path.exists(f'{directory}/modout'):
        # did run already
        return 'alreadyran', '', ''

    out, err, reason = STARS.run(timeout=timeout,
                                 cwd=directory)

    # check last model -- out file and sneplot
    # check sneplot against out
    # classify each model using perl script

    if reason == 'finished' and do_he_flash == True:
        outfile = STARS.get_last_converged_model(f'{directory}/out',
                                                 as_obj=True)

        classification = kaitiaki.classify.go(f'{directory}/out',
                                              f'{directory}/plot')

        if classification == 8.0:
            try:
                STARS.evolve_through_helium_flash(timeout=timeout,
                                                  filedir=directory)
            except Exception as e:
                reason = ('heflashnotpossible', e)

    return out, err, reason


class STARSController:
    """Represents an instance of the STARS code."""
    def __init__(self, verbose_output: bool = True, run_bs: str = "."):
        """Represents an instance of the STARS code."""
        self._verbose_output = verbose_output
        self._output_dir = "."
        self._params = dict()
        self._run_bs_location = run_bs
        self._datafile = "data"

    def blit(self, ZS='z020', directory='.'):
        """Blits the directory given by creating the COTables and data file.

        Loads the default data file and the default COtable.

        Args:
            ZS: the metallicity (BPASS-formatted) to load opacity tables for.
            directory: The directory to blit.

        Returns:
            None

        Notes:
            - Also sets the parameters ZS and CH in the data file. CH is set to 0.75-2.5*ZS.
            - Creates a new file *data* in the directory specified.
            - Creates a new file *COtables* in the directory specified

        """
        # Copy in data, COtable
        data = self.fetch_datafile()
        COtables = kaitiaki.load_file(f'COtables/COtables_{ZS}')

        with open(f'{directory}/data', 'w') as f:
            f.write(data)
        with open(f'{directory}/COtables', 'w') as f:
            f.write(COtables)

        if 'em' not in ZS:
            Z = '0.' + ZS[1:]
        else:
            Z = '1e-' + ZS[-1]

        Z = float(Z)

        self.configure_parameters({'ZS': Z, 'CH': 0.75-2.5*Z})

    def set_period(self, period,
                         directory='.',
                         boost_max_nmodels=True,
                         forcibly_do_both=False):
        """Sets the period in the modin file.

        """
        modin_location = directory + "/modin"
        modins = [modin_location]

        in_secondary_mode = ('imode' in self._params.keys() and self._params['imode'] == 2)

        if in_secondary_mode or forcibly_do_both:
            # We are in binary evolution mode
            modins.append(f"{modin_location}2")

        for modin in modins:
            with open(modin, 'r+') as f:
                data = f.readlines()

                before_period = data[0][:46]
                after_period = data[0][60:]
                period = '{:12.6E}'.format(float(period))
                string = before_period + '  ' + period + after_period

                data[0] = string

                if boost_max_nmodels:
                    data[0] = data[0][:94] + " 99999      0" + data[0][107:]

                f.seek(0)
                f.writelines(data)
                f.truncate()

    def update_datafile(self, new_location):
        self._datafile = new_location

    def update_run_bs(self, loc):
        self._run_bs_location = loc

    def fetch_datafile(self):
        return kaitiaki.load_file(f"data.bak")

    def load_default_modin(self,
                           directory='.',
                           as_secondary=False,
                           Z='z020',
                           set_nmodels_to=99999):

        dest_file = 'modin'

        if as_secondary:
            dest_file = 'modin2'

        with open(f'{directory}/{dest_file}', 'w') as file:
            modin = kaitiaki.load_file(f'modins/modin.bak.{Z}')
            modin = modin.split("\n")

            if set_nmodels_to != None:
                nmods = str(set_nmodels_to).rjust(5, '0')
                modin[0] = modin[0][:94] + f" {nmods}" + modin[0][100:]

            file.write("\n".join(modin))

    # TODO: Fix (dir isn't changed correctly and this will also
    # probably not actually put the file anywhere
    def generate_datafile(self, loc):
        dfile = self.fetch_datafile()
        wd = os.getcwd()

        for d in loc.split("/"):
            os.chdir(d)

        with open(loc, 'w') as f:
            f.write(dfile)

    def output(self, msgtype, message):
        """Outputs a message to stdout, if verbose_output is on. Does
        nothing otherwise.

        Args:
            msgtype (str): The message type (warning, error, info, status) to output
            message (str): The message to output
        """
        if self._verbose_output:
            kaitiaki.debug(msgtype, str(message))

    def configure_parameters(self, params):
        """Sets up the parameters. Doesn't commit them yet.

        Arguments:
            params {dict} -- A dictionary containing key: value pairs
            of what to write to the datafile.
        """
        for k, v in params.items():
            self._params[k.lower()] = v

    def _write_orbital_equations(self, block, dfile_name=''):
        if dfile_name == '':
            dfile = self._datafile

        with open(dfile_name, 'r+') as f:
            dfile = f.readlines()

            f.seek(0)

            blk = block.split("\n")
            for i in range(4):
                dfile[4+i] = blk[i] + "\n"

            f.writelines(dfile)
            f.truncate()

    def setup_binary_evolution(self, dfile='data'):
        """
        Modifies data to allow for binary evolution. Sets the following
        parameters:
            ID block for binaries
            IMODE   - To 2 (binaries)
            IML1    - To 5 (custom) TODO: Check what prescription is
            IML2    - To 5 (custom)
            RML     - 0 (off)
            ITH     - 1 (on)
            IX      - 1 (on)
            IY      - 1 (on)
            IZ      - 1 (on)
            ISTART  - 1 (reset age, nmod, dt)
        """
        binary_block = """ 14 14  0  9  1102  0  0  0 99
  1  2  4 16 17 19 13 14 29  5  3  9 10 11 12 15 20 18 24 25 26 27 30  8  7  6 23 22 21  0
  7  8  9 10 11 12 14 22 23 24 25 26 27 29  4  2  1  3  5  6 19 17 16 13 28 18 20 21  0  0
  4  5  6  7  8  9 10 19 20 21 22 23 24 25  2  3  1 17 18 16  2  3  1  4  5 17 18 16 20  0"""

        self._write_orbital_equations(binary_block, dfile)

        params = {
            'imode': 2
        }

        self.configure_parameters(params)

    # TODO: I have to make this fix the period if you set a single star
    # after a binary run
    def setup_single_evolution(self, dfile='data'):
        """
        Modifies data to allow for single star evolution. Sets the following
        parameters:
            ID block for single stars
            IMODE   - To 1 (single stars)
            IML1    - To 5 (custom) TODO: Check what prescription is
            IML2    - To 0 (off - shouldn't matter though)
            RML     - 0 (off)
            ITH     - 1 (on)
            IX      - 1 (on)
            IY      - 1 (on)
            IZ      - 1 (on)
            ISTART  - 1 (reset age, nmod, dt)
        """
        single_block = """  6  7  0  3  0 76  0  0  0 99
  1  2  4  5  3  9 10 11 12 15  8  7  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  7  8  9 10 11 12 14  4  2  1  3  5  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  4  5  6  7  8  9 10  2  3  1  2  3  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"""

        self._write_orbital_equations(single_block, dfile)

        params = {
            'imode': 1
        }

        self.configure_parameters(params)

    def setup_zams_inflation(self, mass):
        """
        Sets up the data and modin files for ZAMS inflation.

        Arguments:
            mass {float} -- the ZAMS mass to set.

        Notes:
            This *does not* actually inflate the star -- you still must call
            STARS.run() after this.

        """
        params = {
            'IML1': 9,
            'IML2': 0,
            'RML': float(mass),
            'IMODE': 1,
            'IX': 0,
            'IY': 0,
            'IZ': 0,
            'ITH': 0
        }

        self.setup_single_evolution(dfile=self._datafile)
        self.configure_parameters(params)

    def setup_evolution(self, mode='single', iml_option=5):
        match mode:
            case 'single':
                self.setup_single_evolution()
                IML1 = iml_option
                IML2 = 0
            case 'binary':
                self.setup_binary_evolution()
                IML1 = iml_option
                IML2 = iml_option
            case _:
                error_message = (f"mode must be 'single' or 'binary', "
                                 f"not '{mode}'.")

                raise ValueError(error_message)

        self.configure_parameters({
            'IX': 1,
            'IY': 1,
            'IZ': 1,
            'ITH': 1,
            'ISTART': 1,
            'IML1': IML1,
            'IML2': IML2
            })

    def relax_model(self):
        self.output('status', "Relaxing current modin...")

        data = kaitiaki.file.data(self._datafile)

        self.setup_single_evolution()

        old_params = data.get(['IML1', 'IML2', 'ITH', 'IX', 'IY', 'IZ', 'NWRT1'])

        self.configure_parameters({
            'IML1': 0,
            'IML2': 0,
            'ITH': 0,
            'IX': 0,
            'IY': 0,
            'IZ': 0,
            'NWRT1': 1,
            'ISTART': 1
            })

        self.run()

        self.configure_parameters(old_params)

        self.output('status', 'Model relaxed. It is in modout -- use modout_to_modin() to load it for the next run.')

    def set_output_directory(self, directory):
        """Sets up the output directory.

        If the directory doesn't exist, this will create it.

        Arguments:
            directory {string} -- The directory to output files to.
        """
        self._output_dir = directory

        if not os.path.exists(self._output_dir):
            self.terminal_command(f'mkdir {self._output_dir}')

    def terminal_command(self, command, timeout=5*60, cwd=None, warn=True):
        """
        Included for backwards compatibility.

        You should use kaitiaki.terminal.execute() instead
        """
        return kaitiaki.terminal.execute(command, timeout, cwd, warn)

    def commit_parameters(self):
        with kaitiaki.file.data(self._datafile) as dfile:
            dfile.backup_if_not_exists()

            for param, val in self._params.items():
                dfile.set(param, val)

    def get_binary_mode(self, path):
        with kaitiaki.file.data(f"{path}/data") as data:
            IS_BINARY = bool(int(data.get('imode'))-1)
        return IS_BINARY

    def live_hr_diagram(self, the_dir, compute_thread, fps=60):
        path = os.path.realpath(the_dir)
        IS_BINARY = self.get_binary_mode(path)

        old_params = plt.rcParams
        params = {'text.color': '(0.2,0.2,0.2)',
                  'axes.edgecolor': '(0.1,0.1,0.1)',
                  'axes.labelcolor': '(0.1,0.1,0.1)',
                  'axes.linewidth': '1',
                  'axes.spines.top': 'True',
                  'axes.spines.right': 'True',
                  'axes.spines.bottom': 'True',
                  'axes.spines.left': 'True',
                  'axes.grid': 'True',
                  'grid.alpha': '0.7',
                  'grid.linestyle': '--',
                  'grid.linewidth': '0.6',
                  'axes.prop_cycle': "cycler('color', ['EE6123', '3185FC', 'F49CBB', '694A38', 'FDCA40', '26C485', 'A3E7FC', '32908F', '553A41', '2F0601', '50514F', 'FFE066', '247BA0', '70C1B3'])",
                  'lines.linewidth': '1.5',
                  'lines.markeredgewidth': '0.0',
                  'scatter.marker': '.',
                  'xtick.color': '(0.2,0.2,0.2)',
                  'xtick.direction': 'in',
                  'ytick.color': '(0.2,0.2,0.2)',
                  'ytick.direction': 'in',
                  'xtick.top': 'True',
                  'xtick.bottom': 'True',
                  'ytick.left': 'True',
                  'ytick.right': 'True',
                  'xtick.minor.visible': 'True',
                  'ytick.minor.visible': 'True',
                  'font.family': 'serif',
                  'font.serif': 'CMU Serif',
                  'font.sans-serif': 'CMU Sans Serif',
                  'mathtext.fontset': 'cm',
                  'mathtext.default': 'sf',
                  'figure.labelsize': 'large',
                  'figure.titlesize': 'large',
                  'font.size': 12.0,
                  'xtick.labelsize': 'medium',
                  'ytick.labelsize': 'medium'}

        plt.rcParams.update(params)

        def get_quantity(dataArray, column_number):
            returnable = []

            for eachLine in dataArray:
                if len(eachLine) > 1:
                    line_parts = eachLine.split()

                    quantity = line_parts[column_number].strip()

                    try:
                        if column_number == 0:
                            quantity = int(quantity)
                        else:
                            quantity = float(quantity)

                        returnable.append(quantity)
                    except ValueError:
                        returnable.append(np.nan)

            return returnable

        do_limits = False

        for j, fname in enumerate(['plot', 'plot2']):
            try:
                full_fname = f"{path}/{fname}"
                with open(full_fname, "r") as file:
                    file_data = file.read()
            except FileNotFoundError:
                continue

            dataArray = file_data.split('\n')

        hrd_figure = plt.figure(figsize=(12, 12))

        ax = hrd_figure.add_subplot(111)

        hr_axis = ax  # HR Diagram

        hr_axis.grid(True)

        def format_age(age):
            returned_age = None
            if age < 1e3:
                returned_age = f"{age:.3f} years"
            if age <= 1:
                returned_age = f"{age*365.25:.3f} days"
            if age <= 1/365.25:
                returned_age = f"{age*365*24:.3f} hours"
            if age <= 1/(365.25*24):
                returned_age = f"{age*365*24*24:.3f} minutes"
            if age <= 1/(365.25*24*24):
                returned_age = f"{age*365*24*24*60:.3f} seconds"
            if age >= 1e3:
                returned_age = f"{age/1e3:.3f} kyr"
            if age >= 1e6:
                returned_age = f"{age/1e6:.3f} Myr"
            if age >= 1e9:
                returned_age = f"{age/1e9:.3f} Gyr"
            if age >= 1e12:
                returned_age = f"{age/1e12:.3f} Tyr"
            if age < 0:
                returned_age = f"-{format_age(np.abs(age))}"

            return returned_age

        def animate(i):
            nonlocal IS_BINARY, compute_thread, old_params, hrd_figure

            if not compute_thread.is_alive():
                # Computation finished -- kill this process!
                plt.rcParams.update(old_params)
                plt.close(hrd_figure)
                raise StopIteration

            ax.clear()

            def RL(q, a):
                return a * 0.49 * q**(2/3) / (0.6*q**(2/3) + np.log(1+q**(1/3)))

            def RL2(q, a):
                # https://ui.adsabs.harvard.edu/#abs/2016A&A...588A..50M/abstract
                return (1+0.299*np.arctan(1.849*q**0.397)*q**0.520) * RL(q, a)

            cs = ['purple', 'skyblue']
            labels = ['Primary', 'Secondary']

            models = [0]  # to block the legend on first run.

            for j, fname in enumerate(['plot', 'plot2']):
                try:
                    with open(f"{path}/{fname}", "r") as file:
                        file_data = file.read()
                except FileNotFoundError:
                    continue

                dataArray = file_data.split('\n')

                temperature = get_quantity(dataArray, 3)
                luminosity = get_quantity(dataArray, 4)
                models = get_quantity(dataArray, 0)
                ages = get_quantity(dataArray, 1)
                age = format_age(ages[-1])

                if len(temperature) < 10:
                    # Probably a new run
                    IS_BINARY = self.get_binary_mode(path)

                if IS_BINARY:
                    suffix = "Binary"
                else:
                    suffix = "Single"

                if len(ages) > 2:
                    dt = format_age(ages[-1] - ages[-2])
                else:
                    dt = format_age(ages[-1])

                title = f"""Evolution mode: {suffix}
Model Number {models[-1]}
Current Age: {age}
dt: {dt}"""

                hr_axis.text(0.01, 0.99, title,
                             ha='left',
                             va='top',
                             transform=hr_axis.transAxes,
                             fontsize=14)

                hr_axis.plot(temperature, luminosity, c=cs[j], label=labels[j])

                hr_axis.scatter(temperature[::1000], luminosity[::1000],
                                marker='o',
                                edgecolors='orange',
                                facecolors='none')

                hr_axis.set_xlabel(r"$\log~\mathrm{T}_\text{eff}~/~\text{K}$")
                hr_axis.set_ylabel(r"$\log~\mathrm{L}~/~\text{L}_\odot$")

            if models[-1] > 0:
                hr_axis.legend(frameon=False)

            if do_limits:
                hr_axis.set_xlim(4.7, 3.5)
                hr_axis.set_ylim(3.5, 6.0)
            else:
                if hr_axis.get_xlim()[0] < hr_axis.get_xlim()[1]:
                    hr_axis.axes.invert_xaxis()

        try:
            ani = animation.FuncAnimation(hrd_figure, animate,
                                          repeat=False,
                                          interval=1000/fps)
            plt.show()
        except StopIteration:
            return

    def run(self,
            timeout=5*60,
            cwd=None,
            warn=True,
            time_me=False,
            with_live_HR=False,
            HRD_fps=60):
        """Runs ./run_bs

        Returns:
            tuple -- The output from run() above.
        """
        wd = os.getcwd()
        dirname = wd

        if cwd is not None:
            os.chdir(cwd)
            dirname = cwd

        self.commit_parameters()

        if time_me:
            time_start = time.time_ns()

        cmd = f'{self._run_bs_location}/run_bs'

        if with_live_HR:
            ###
            # Known bug:
            # This will cause the code to segfault.
            ###

            # Main thread owns the HRD, background thread owns stars.

            compute = threading.Thread(target=self.terminal_command,
                                       name="AotearoaSTARS",
                                       args=(cmd, timeout, dirname, warn))
            compute.start()
            try:
                self.live_hr_diagram(dirname, compute, fps=HRD_fps)
            except AttributeError:
                pass  # It throws an attributeerror when it stops...
                      # for some reason.
        else:
            out, err, reason = self.terminal_command(cmd,
                                                     timeout=timeout,
                                                     warn=warn)

            if time_me:
                time_end = time.time_ns()

                delta_time = time_end - time_start

            os.chdir(wd)

            if time_me:
                return out, err, reason, delta_time
            else:
                return out, err, reason

    def get_last_converged_model(self, file, as_obj=False):
        from file_read_backwards import FileReadBackwards
        from queue import Queue

        q = []

        i = 0

        with FileReadBackwards(file, encoding="utf-8") as frb:
            for line in frb:
                while len(q) > 8:
                    q.pop(0)
                q.append(line)
                if 'dt/age/MH/MHe' in line.strip():
                    break

        modelblock = []

        for line in reversed(q):
            modelblock.append(line)

        if not as_obj:
            return modelblock

        return kaitiaki.file.outfile.ModelSummary(modelblock)

    def get_outfile_models(self, file: str):
        models = []
        try:
            # this is a file path
            with open(file, 'r') as out:
                contents = out.readlines()
                outlen = len(contents)
                for i in range(outlen):
                    if 'dt/age/MH/MHe' in contents[i]:
                        # we got one
                        modelnum = contents[i].strip().split()[0]
                        MHe = contents[i+2].strip().split()[0]
                        mass = contents[i+1].strip().split()[0]
                        model = contents[i:i+8]
                        models.append((modelnum, mass, MHe, model))
                    i += 8
            return models
        except OSError:
            # this is the contents of a file
            contents = file.split("\n")
            outlen = len(contents)
            for i in range(outlen):
                if 'dt/age/MH/MHe' in contents[i]:
                    # we got one
                    modelnum = contents[i].strip().split()[0]
                    MHe = contents[i+2].strip().split()[0]
                    mass = contents[i+1].strip().split()[0]
                    model = contents[i:i+8]
                    models.append((modelnum, mass, MHe, model))
                i += 8
            return models

    def evolve_async(self,
                     masses,
                     timeout: int = 20*60,
                     evolution_dir: str = "",
                     attempt_he_flash: bool = True,
                     data_params: dict = dict(),
                     reserve_core: bool = True,
                     logged_masses: bool = False,
                     ZAMS_files_location: str = "",
                     metallicity: str = "z020"):
        # OLD_DFILE_LOC = self._datafile

        if len(data_params.keys()) == 0:
            # use default?
            params = dict()
            self.output('warning', "data_params not set!")
        else:
            # use these ones?
            params = data_params

        cpu = _allocate_cores(reserve_core)
        pool = mp.Pool(cpu)

        outputs = []
        outs = []

        with tqdm(total=len(masses)) as pbar:
            for mass in masses:
                directory = f"{evolution_dir}/{metallicity}/m{mass}"

                # ok first thing, copy data down here:
                data = kaitiaki.load_file(f"data.bak")

                with open(f"{directory}/data", "w") as f:
                    f.write(data)
                    # self._datafile = f"{directory}/data"

                # we also need to modify the data file.
                with kaitiaki.file.data(f"{directory}/data") as dfile:
                    for param, val in params.items():
                        dfile.set(param, val)
                    self._params = dict()

                # now copy in modin
                self.modout_to_modin(f"{ZAMS_files_location}/m{mass}/modout",
                                     f"{directory}/modin")

                # work out how many parent directories we have to go up
                # to find run_bs.
                num_run_bs = len(directory.rstrip("/").split("/"))
                run_bs = ("../" * num_run_bs).rstrip("/")

                if logged_masses:
                    mass = 10**mass

                def update(*a): pbar.update(1)
                out = pool.apply_async(_worker_evolve,
                                       args=(directory,
                                             timeout,
                                             mass,
                                             True,
                                             run_bs,
                                             self),
                                       callback=update)

                outputs.append((mass, out))

            for r in outputs:
                outs.append((r[0], r[1].get()))

        # self._datafile = OLD_DFILE_LOC

        return outs

    def evolve_through_helium_flash(self, timeout=30*60, filedir="."):
        modelblock = self.get_last_converged_model(f'{filedir}/out')

        required_he_core_mass = modelblock[2].strip().split()[0]
        target_mass = modelblock[1].strip().split()[0]

        out = kaitiaki.load_file(f"pseudo_evolution/out")
        models = self.get_outfile_models(out)

        def near(x, y, threshold=0.01):
            x = float(x)
            y = float(y)

            return x in kaitiaki.helpers.Range(y-threshold*y, y+threshold*y)

        try:
            he_mass = required_he_core_mass

# our pre-generated 3Msun model starts at NMOD=5000.
            mod_nums = {int(mod[0]) for mod in models if near(mod[2], he_mass)}
            new_model_number = min(mod_nums)# - 5000

            self.output("status", f"He Core Mass required: {he_mass}")
        except ValueError:
            self.output("error", (f"No matching model found for He "
                                  f"core mass {required_he_core_mass}. "
                                  f"The requested model mass was "
                                  f"{target_mass} Msun."))

            raise ValueError((f"No matching model found for He "
                              f"core mass {required_he_core_mass}. "
                              f"The requested model mass was "
                              f"{target_mass} Msun."))

        self.output("status", "Loading the shipped modout file. This is 1.5GB so it may take a second!")

        # Extract this model, write it to modin.
        # modout = _load_shipped_file(f"../backup_data/pseudo_evolution/modout")
        # self.output("status", "File loaded.")

        from backup_data import pseudo_evolution as resources
        modout = pkg_resources.open_text(resources, 'modout')

        self.output("status", "Determining the relevant model...")

        with open(f"{filedir}/modin", "w") as modin_f:
            nmesh = int(modout.readline()[88:88+6].strip())*2+1
            self.output('status', f'nmesh is: {nmesh}')

            nlines = (nmesh-1)//2
            self.output('status', f'nlines is: {nlines}')

            i = int(nmesh*new_model_number)
            self.output('status', f'i is: {i}')

            model = []

            for j, line in enumerate(modout):
                if i-2 == j:
                    # this is the start of the block.
                    # Now iterate over the rest of the model file
                    # which is the next NMESH lines.
                    for k in range(nmesh):
                        line = modout.readline()  # grab the line
                        model.append(line)  # store it
                    break  # ok, we only need this model, so stop here.

            self.output('status', f'model found. writing...')

            modin_f.writelines(model)

        modout.close()

        new_nmod = new_model_number

        self.output('status', (f"Inferred required model number from "
                               f"pre-generated 3Msun star: {new_nmod}"))

        new_model_number = int(modelblock[0].strip().split()[0])# - 5000

        self.output('status', (f"New model number from the out file you "
                               f"just ran: {new_model_number} (should be "
                               "close to last printout)"))

        self.output('info', f'Prepared the new modin')

        with kaitiaki.file.data(self._datafile) as dfile:
            dfile.backup_if_not_exists()

        for file in ['modout', 'out', 'plot']:
            self.terminal_command((f"cp {filedir}/{file} "
                                   f"{filedir}/{file}.prehf"))

        params = {
            'IMODE': 1,          # Single-star mode
            'IML1': 9,           # RE-ML: target a specific mass.
            'IML2': 0,           # No ML for the secondary
                                 # (because there isn't one.)
            'RML': target_mass,  # This is the target ZAMS mass in Msun
            'RMG': 0,            # Outdated, otherwise the code doesn't
                                 # iterate past one loop
            'ITH': 0,            # Thermal expansion off
            'IX': 0,             # Hydrogen burning off
            'IY': 0,             # Helium burning off
            'IZ': 0,             # Metal burning off
            'NCH': 1,
            'ISTART': 1,         # Don't reset age, dt, nmod,
            'NNMOD': new_model_number
        }

        self.configure_parameters(params)

        stdout, stderr, status = self.run(cwd=filedir,
                                          timeout=timeout)

        with open(f'{filedir}/stdout', 'w') as fout:
            with open(f'{filedir}/stderr', 'w') as ferr:
                fout.write(stdout)
                ferr.write(stderr)

        self.output('info', (f'Finished the inflationary stage of the '
                             f'reference model. Reason: {status}'))

        self.output('info', (f'Evolved the reference model to the target '
                             f'mass of {target_mass}'))

        if status == 'finished':
            self.modout_to_modin(modout_location=f"{filedir}/modout",
                                 modin_location=f"{filedir}/modin")

            # self.terminal_command("cp pseudoev/modout.3 modin")

            params = {
                'IMODE': 1,  # Single-star mode
                'IML1': 5,
                'IML2': 0,   # No ML for the secondary (because there is none.)
                'RML': 0,    # This is the target ZAMS mass in Msun
                'RMG': 0,    # Outdated, otherwise the code doesn't
                             # iterate past one loop
                'ITH': 1,
                'IX': 1,
                'IY': 1,     # Helium burning off
                'IZ': 1,     # Metal burning off
                'ISTART': 1,
                'NCH': 2,
                'NNMOD': new_model_number
            }

            self.configure_parameters(params)

            stdout, stderr, status = self.run(cwd=filedir, timeout=timeout)

            self.output('info', f'Completed evolution. Now moving files.')

            if status == 'timeout':
                self.output('warning', ("run_bs timed out on the "
                                           "second evolution run. Output "
                                           "has been saved but consider "
                                           "rerunning with a longer "
                                           "timeout."))

            if status in ['timeout', 'finished']:
                self.terminal_command((f"cp {filedir}/modout "
                                       f"{filedir}/modout.posthf"))

                self.terminal_command((f"cp {filedir}/plot "
                                       f"{filedir}/plot.posthf"))

                self.terminal_command((f"cp {filedir}/out "
                                       f"{filedir}/out.posthf"))
            else:
                self.output('error', f"An error occurred in the second evolution stage. Reason: {status}")
                self.output('error', stderr)
        else:
            self.output('error', f"An error occurred in the first evolution stage. Reason: {status}")
            self.output('error', stderr)

    def modout_to_modin(self,
                        modout_location="modout",
                        modin_location="modin",
                        n_models_ago=0):
        """Moves the last model in modout to modin

        Reads modout to determine the number of lines to move (columns 88-94 of the first line).

        Args:
            modout_location (str): The modout file to load (default: `"modout"`)
            modin_location (str): The modin file to write to (default: `"modin"`)
            n_models_ago (number): How many models ago to load (e.g. n_models_ago=3 loads 3 models back, n_models_ago=-1 loads the first model in the file, n_models_ago=0 loads the last model in the file.) (default: `0`)
    """
        with open(modout_location, 'r') as f:
            nmesh = int(f.readline()[88:94])

        lines = nmesh * 2 + 1

        with open(modout_location, 'r') as f:
            output = f.readlines()

        N = len(output)

        if n_models_ago >= 0:
            output = output[N-(n_models_ago+1)*lines:N-(n_models_ago)*lines]
        else:
            n_models_ago = np.abs(n_models_ago)
            output = output[(n_models_ago-1)*lines:(n_models_ago)*lines]

        with open(modin_location, 'w') as f:
            f.writelines(output)

    def run_default_evolution(self, zams_mass, NM2=199):
        """Performs one run from pre-ZAMS until the end of evolution.
        Assumes NORMAL file structure.

        Not really useful for scientific applications -- just here
        for educational purposes on how STARS runs.

        Arguments:
            zams_mass {float} -- The ZAMS mass we should target.
        """
        self.output('status', 'Preparing for ZAMS evolution')
        self.terminal_command("rm modin")
        self.terminal_command("rm data")

        self.blit()
        self.load_default_modin()

        pre_zams_params = {'IML1': 9,
                           'RML': float(zams_mass),
                           'ITH': 0,
                           'IX': 0,
                           'IY': 0,
                           'IZ': 0,
                           'NM2': NM2,
                           'NCH': 3,
                           'ISTART': 0}

        self.output('status', 'Configuring Parameters')
        self.configure_parameters(pre_zams_params)
        self.output('status', 'Beginning ZAMS evolution')
        output, error, status = self.run(timeout=20*60)

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
                           'ISTART': 1,
                           'NM2': NM2,
                           'NNMOD': 0}

        self.output('status', 'Configuring Parameters')
        self.configure_parameters(main_evo_params)

        self.output('status', 'Beginning post-ZAMS evolution')

        output, error, status = self.run(timeout=20*60)

        if self._verbose_output:
            self.output("status", output)
            self.output("error", error)

        self.output('status', 'Moving output files')

        self.terminal_command(f'mv plot {out}/plot.model')
        self.terminal_command(f'mv out {out}/out.model')
        self.terminal_command(f'mv modout {out}/modout.model')
