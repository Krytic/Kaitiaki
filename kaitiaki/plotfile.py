from copy import deepcopy
from decimal import Decimal
from os import path

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

class Plotfile:
    def __init__(self, file, attempt_loading_mplstylefile=False):
        self._data, status = self.parse_plotfile(file)
        self._segment_points = []

        if status == 'skipped':
            raise OSError("Requested file does not exist")

        # This is a staggered attempt to load an inbuilt matplotlib
        # style file.
        if attempt_loading_mplstylefile:
            try:
                # See if the krytic stylefile is available...
                plt.style.use("krytic")
            except: # it isn't...
                try:
                    # what if it's in the current directory?
                    plt.style.use("krytic.mplstyle")
                except: # it isn't...
                    # Screw it, we'll just ignore it (it's a style file,
                    # it doesn't matter for actual plotting, just
                    # visualisation)
                    pass

    def __add__(self, other):
        if isinstance(other, Plotfile):
            self._segment_points.append(len(self._data))
            self._data = self._data.append(other._data)

            return self
        else:
            raise TypeError("One of these items isn't a Plotfile object.")

    def access(self):
        return self._data

    def hr_diagram(self, **kwargs):
        self.plot('logT', 'logL', **kwargs)
        plt.gca().invert_xaxis()

    def plot(self, x_axis, y_axis, **kwargs):
        """Plots the parameter given by x_axis against y_axis.

        Arguments:
            x_axis {str} -- The x-axis to plot. Must be a key of self._data
            y_axis {str} -- The y-axis to plot. Must be a key of self._data
            **kwargs {variable} -- Any keyword arguments to be passed to plt.plot.

        Returns:
            obj {list} -- The list of Line2d objects plotted by plt.plot.
        """
        x_arr = self._data[x_axis].to_numpy()
        y_arr = self._data[y_axis].to_numpy()

        # Check if we have to stitch together multiple files
        if len(self._segment_points) > 1:
            for i in range(len(self._segment_points)):
                si = self._segment_points[i]
                if i == 0:
                    x = x_arr[0:si]
                    y = y_arr[0:si]
                elif i < len(self._segment_points):
                    li = self._segment_points[i-1]
                    x = x_arr[li:si]
                    y = y_arr[li:si]
                else:
                    x = x_arr[si:]
                    y = y_arr[si:]

                obj = plt.plot(x, y, **kwargs)
        else:
            x = x_arr
            y = y_arr
            obj = plt.plot(x, y, **kwargs)

        return obj

    def parse_plotfile(self, fname):
        # https://youtu.be/4THFRpw68oQ?t=34
        c = ['N', 'age', 'logR', 'logT', 'logL', 'M', 'MH', 'MHe', 'LH',
             'LHe', 'LC', 'Mconv_1', 'Mconv_2', 'Mconv_3', 'Mconv_4',
             'Mconv_5', 'Mconv_6', 'Mconv_7', 'Mconv_8', 'Mconv_9',
             'Mconv_10', 'Mconv_11', 'Mconv_12', 'MHmax', 'MHemax',
             'logK', 'dt', 'XHs', 'XHes', 'XCs', 'XNs', 'XOs', 'X3Hes',
             'R/RL', 'J1', 'PBin', 'rsep', 'Mtot', 'Jorb', 'J1+J2', 'Jtot',
             'worb', 'w1', 'I1', 'Iorb', 'Mdot', 'Mshell_1', 'Mshell_2',
             'Mshell_3', 'Mshell_4', 'Mshell_5', 'Mshell_6', 'Mshell_7',
             'Mshell_8', 'Mshell_9', 'Mshell_10', 'Mshell_11', 'Mshell_12',
             'Mth_1', 'Mth_2', 'Mth_3', 'Mth_4', 'Mth_5', 'Mth_6', 'Mth_7',
             'Mth_8', 'Mth_9', 'Mth_10', 'Mth_11', 'Mth_12', 'Mconv_env',
             'Rconv_env', 'logrho', 'logTc']

        if path.exists(f'{fname}'):
            # Following is a python implementation of the following
            # FORTRAN 77 format statement. We must encode this manually.
            # Note that Pandas does have an infer_nrows option for
            # pd.read_fwf, but the inference can be a little strange
            # sometimes, so I prefer to manually define it (since we
            # know the widths a priori anyway from STARS). Note that
            # the P statement adjusts the scaling factor and does not
            # contribute to the output short of determining precisely
            # where the decimal point sits. For this purpose, however,
            # that is irrelevant.
            # I6,1P,E16.9,0P,24F10.5,1P,3E13.6,18(1X,E12.5),0P,52F9.5
            spec = ([6,16]                               # I6, E16.9
                 + [10 for _ in range(24)]               # 24F10.5
                 + [13, 13, 13]                          # 3E13.6
                 + [i for _ in range(18) for i in (1,12)]#18(1X,E12.5)
                 + [9 for _ in range(52)])               # 52F9.5
            # The spec extends out to ~100 columns to future proof it
            # I think, so we have to truncate it here to the length of
            # what we know is in the file.
            spec = spec[:len(c)]

            # That giant comment aside, though, something is broken in
            # that spec because some values are not being read properly
            # in some test files (e.g., reads NaN instead of 1e+34).

            # Temporary fix: set infer_nrows to 99999, the highest
            # theoretical length of plot.

            df = pd.read_fwf(fname,
                             names=c,
                             infer_nrows=99999)

            status = 'loaded'
        else:
            df = None
            status = 'skipped'

        return df, status

class ClusteredPlotfile:
    def __init__(self, directory, indexes, file_names=['zams', 'model']):
        self._pfiles = dict()


        for idx in indexes:
            n = 0

            for file in file_names:
                F = Plotfile(f'{directory}/plot.{file}.{idx}')

                if n == 0:
                    pfile = F
                else:
                    pfile = pfile + F

                n += 1

            self._pfiles[idx] = pfile

        self._vmin = np.amin(indexes)-0.5
        self._vmax = np.amax(indexes)+0.5

    def plot(self, x_axis, y_axis, cbar_label, cmap='jet', alpha=1):
        fig = plt.figure()

        norm = matplotlib.colors.Normalize(
            vmin=self._vmin,
            vmax=self._vmax)

        # choose a colormap
        c_m = plt.get_cmap(cmap)

        # create a ScalarMappable and initialize a data structure
        s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
        s_m.set_array([])

        for idx, pfile in self._pfiles.items():
            pfile.plot(x_axis, y_axis,
                       color=s_m.to_rgba(idx),
                       alpha=alpha)

        fig.subplots_adjust(right=0.8)
        # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(s_m)#, ax=fig.gca())

        cbar.set_label(cbar_label)

        return fig