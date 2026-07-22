import kaitiaki
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# should move to helpers.py?
def _find_runs(x):
    """Find runs of consecutive items in an array."""

    # ensure array
    x = np.asanyarray(x)
    if x.ndim != 1:
        raise ValueError('only 1D array supported')
    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_starts, run_lengths


class Star:
    def __init__(self, files):
        self.__plot = kaitiaki.file.plot(files['plot'])
        self.__out = kaitiaki.file.out(files['out'])

    def set_out(self, new_out):
        self.__out = new_out

    def get_out(self):
        return self.__out

    def plot(self, *args, **kwargs):
        return self.__plot.plot(*args, **kwargs)

    def get(self, *args, **kwargs):
        return self.__plot.get(*args, **kwargs)

    def plotfile(self):
        return self.__plot

    def core_exhaustion(self, element, threshold=1e-5):
        ind, abundance = self.core_abundance(element)

        try:
            idx = ind[np.where(abundance < threshold)[0][0]]
            return idx, True
        except IndexError:
            return 0, False

    def core_abundance(self, element):
        return self.abundance_at_meshpoint(element, 0)

    def surface_abundance(self, element):
        return self.abundance_at_meshpoint(element, self.__out.get('NM2')-1)

    def abundance_at_meshpoint(self, element, meshpoint):
        ind, NM, res = self.__out.build_network(element)

        return ind, res[:, meshpoint]


class Binary:
    def __init__(self, primary_data, secondary_data):
        if isinstance(primary_data, Star):
            self.__primary = primary_data
            self.__secondary = secondary_data
        else:
            self.__primary = Star(primary_data)
            self.__secondary = Star(secondary_data)

    def plot(self, *args, **kwargs):
        fig, axes = plt.subplots(1, 2)

        self.__primary.plot(*args, **kwargs, ax=axes[0])
        self.__secondary.plot(*args, **kwargs, ax=axes[1])

        return fig, axes

    def get_primary(self):
        return self.__primary

    def get_secondary(self):
        return self.__secondary

    def contact_phases(self):
        donor = kaitiaki.helpers.roche_lobes(self.__primary, self.__secondary)
        accretor = kaitiaki.helpers.roche_lobes(self.__secondary, self.__primary)

        donor_L1 = donor['L1'] <= 10**self.__primary.get('log(R)')
        accretor_L1 = accretor['L1'] <= 10**self.__secondary.get('log(R)')

        in_contact = np.logical_and(donor_L1, accretor_L1).astype(int)

        rv, rs, rl = _find_runs(in_contact)

        phases = []

        for i in range(len(rv)):
            if int(rv[i]) == 1:
                phases.append((rs[i], rs[i]+rl[i]))

        return phases

    def hr_diagram(self, *args, **kwargs):
        fig, axes = plt.subplots(1, 1)

        self.__primary.plotfile().hr_diagram(*args, **kwargs, ax=axes)
        self.__secondary.plotfile().hr_diagram(*args, **kwargs, ax=axes)

        return fig, axes

    def kippenhahn_diagram(self, *args, **kwargs):
        fig, axes = plt.subplots(1, 2)

        self.__primary.plotfile().kippenhahn_diagram(*args,
                                                     **kwargs,
                                                     ax=axes[0])

        self.__secondary.plotfile().kippenhahn_diagram(*args,
                                                       **kwargs,
                                                       ax=axes[1])

        return fig, axes


def load_directory(directory='.', fext=''):
    if fext != '' and not fext.startswith('.'):
        fext = f'.{fext}'

    primary = {
        'out': f'{directory}/out{fext}',
        'plot': f'{directory}/plot{fext}'
    }

    secondary = {
        'out': f'{directory}/out2{fext}',
        'plot': f'{directory}/plot2{fext}'
    }

    return Binary(primary, secondary)
