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
    def __init__(self, file):
        self._data, status = self.parse_plotfile(file)
        self._segment_points = []

        if status == 'skipped':
            raise OSError("Requested file does not exist")

    def __add__(self, other):
        if isinstance(other, Plotfile):
            self._segment_points.append(len(self._data))
            self._data = self._data.append(other._data)

            return self
        else:
            raise TypeError("One of these items isn't a Plotfile object.")

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
          df = pd.read_csv(fname,
                           names=c,
                           sep='\s+')

          status = 'loaded'
        else:
          df = None
          status = 'skipped'

        return df, status

class ClusteredPlotfile:
    def __init__(self, directory, indexes, file_names=['zams', 'model']):
        self._pfiles = dict()

        for idx in indexes:
            zams_file = Plotfile(f'{directory}/plot.zams.{idx}')
            postzams  = Plotfile(f'{directory}/plot.model.{idx}')

            pfile = zams_file + postzams

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