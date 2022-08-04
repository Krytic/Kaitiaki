from copy import deepcopy
from decimal import Decimal
import itertools
from os import path

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

def get_last_line(self, file):
    from file_read_backwards import FileReadBackwards

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

    with FileReadBackwards(file, encoding="utf-8") as frb:
        line = frb.readline()

    spec = ([6,16]                               # I6, E16.9
         + [10 for _ in range(24)]               # 24F10.5
         + [13, 13, 13]                          # 3E13.6
         # should the i in the following LC be a _?
         + [12 for _ in range(18)]# for i in (1,12)]#18(1X,E12.5)
         + [9 for _ in range(52)])               # 52F9.5

    for i, col in enumerate(c):
        width = spec[i]
        # Continue writing parser.

class Plotfile:
    def __init__(self, file,
                       attempt_loading_mplstylefile=False,
                       allow_pad_age=False):

        self._data, status = self.parse_plotfile(file)
        self._segment_points = []

        if status == 'skipped':
            raise OSError("Requested file does not exist")

        self._allow_pad_age = allow_pad_age

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

    def last(self):
        return self._data.iloc[-1:]

    def __add__(self, other):
        if isinstance(other, Plotfile):
            if self._allow_pad_age:
                other.pad_age(self.access()['age'].to_numpy()[-1])

            self._segment_points.append(len(self._data))
            self._data = self._data.append(other._data)


            return self
        else:
            raise TypeError("One of these items isn't a Plotfile object.")

    def access(self):
        return self._data

    def hr_diagram(self, ax=None, **kwargs):
        obj = self.plot('logT', 'logL', ax=ax, **kwargs)
        xlim = obj[0].axes.get_xlim()
        if xlim[0] < xlim[1]:
            obj[0].axes.invert_xaxis()

    def kippenhahn_diagram(self, distinguish_envelopes: bool=False,
                                 legend: bool=True,
                                 x_axis: str='modelnum',
                                 ax=None,
                                 annotate: bool=True,
                                 **kwargs):

        err_msg = "x_axis must be collapsetime, modelnum, or age."

        assert x_axis in ['collapsetime', 'modelnum', 'age'], err_msg

        if ax is None:
            ax = plt.gca()

        if x_axis == 'modelnum':
            X = 'N'
            x_label = "Model Number"
        elif x_axis == 'collapsetime':
            X = 'collapsetime'
            x_label = "Time until collapse [yr]"
        elif x_axis == 'age':
            X = 'age'
            x_label = "Age [yr]"

        self.plot(X, 'MH',
                     c='b',
                     ls='-',
                     label="He core mass",
                     ax=ax) # Helium Core Mass

        self.plot(X, 'MHe',
                     c='k',
                     ls='-',
                     label="CO core mass",
                     ax=ax) # CO Core Mass

        for env in range(1, 13):
            if distinguish_envelopes:
                label = "(semi)conductive envelope"
            else:
                label = "Conductive Envelope"

            lab = None if env < 12 else label

            self.plot(X, f'Mconv_{env}', c='r',
                                         ls='',
                                         marker='.',
                                         label=lab,
                                         absolute=(not distinguish_envelopes),
                                         ax=ax)

        self.plot(X, 'M', c='y', ls='-', label="Total Mass", ax=ax)

        ZAMS = self.access()['M'].to_numpy()[0]

        if annotate:
            ax.set_title(rf"$M_{{\rm ZAMS}}={ZAMS}M_\odot$ star")
            ax.set_xlabel(x_label)
            ax.set_ylabel(r'Mass co-ordinate')

        if legend:
            ax.legend()


    def plot(self, x_axis, y_axis, unlog='neither', age_unit=None, absolute=False, ax=None, **kwargs):
        """Plots the parameter given by x_axis against y_axis.

        Arguments:
            x_axis {str} -- The x-axis to plot. Must be a key of self._data
            y_axis {str} -- The y-axis to plot. Must be a key of self._data
            **kwargs {variable} -- Any keyword arguments to be passed to plt.plot.

        Returns:
            obj {list} -- The list of Line2d objects plotted by plt.plot.
        """
        if x_axis == 'collapsetime':
            age_at_collapse = self.access()['age'].to_numpy()[-1]
            time_until_collapse = age_at_collapse \
                                  - self.access()['age'].to_numpy()

            x_arr = time_until_collapse
        else:
            x_arr = self._data[x_axis].to_numpy()

        y_arr = self._data[y_axis].to_numpy()

        if unlog == 'y':
            y_arr = 10 ** y_arr
        elif unlog == 'x':
            x_arr = 10 ** x_arr
        elif unlog == 'both':
            y_arr = 10 ** y_arr
            x_arr = 10 ** x_arr

        if absolute:
            y_arr = np.abs(y_arr)

        if ax is None:
            ax = plt.gca()
        else:
            ax = ax

        plt.sca(ax)

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
                objs.append(obj)
        else:
            x = x_arr
            y = y_arr
            objs = [plt.plot(x, y, **kwargs)]

        return list(itertools.chain.from_iterable(objs))

    def pad_age(self, by):
        self._data['age'] += by

    def parse_plotfile(self, fname):
        """
        Parses a plotfile.

        NOTE that MHe is the mass of the He-*exhausted* core. So it is
        technically the CO core mass.
        Similarly MH is the H-*exhausted* core! So it is technically the
        He core mass.

        This design decision was made by the maintainers of STARS. We use their
        terminology to be consistent.
        """
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
                 # should the i in the following LC be a _?
                 + [12 for _ in range(18)]# for i in (1,12)]#18(1X,E12.5)
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