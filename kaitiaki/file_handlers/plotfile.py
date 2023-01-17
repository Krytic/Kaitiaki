import kaitiaki

from copy import deepcopy
from decimal import Decimal
import itertools
from os import path
from hoki import load
from hoki import constants as hoki_constants

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

def get_last_line(self, file):
    from file_read_backwards import FileReadBackwards

    c = kaitiaki.constants.PLOT_FILE_COLUMNS

    with FileReadBackwards(file, encoding="utf-8") as frb:
        line = frb.readline()

    spec = ([6,16]                               # I6, E16.9
         + [10 for _ in range(24)]               # 24F10.5
         + [13, 13, 13]                          # 3E13.6
         + [12 for _ in range(18)]               #18(1X,E12.5)
         + [9 for _ in range(52)])               # 52F9.5

    for i, col in enumerate(c):
        width = spec[i]
        # Continue writing parser.

class plot:
    def __init__(self, file='plot',
                       allow_pad_age=True,
                       row='all',
                       dummy_object=False):

        assert row in ['all', 'last'], "row must be all or last"

        self._data, status = self.parse_plotfile(file, row, dummy_object)
        self._segment_points = []

        if status == 'skipped':
            raise OSError("Requested file does not exist")

        self._allow_pad_age = allow_pad_age

    def __len__(self):
        return len(self._data)

    def last(self):
        return self._data.iloc[-1:]

    def zams(self):
        return self._data.iloc[0]

    def __add__(self, other):
        if isinstance(other, Plotfile):
            if self._allow_pad_age:
                other.pad_age(self.get('age').to_numpy()[-1])

            # self._segment_points.append(len(self._data))
            # self._data = self._data.append(other._data)

            new_dataframe = Plotfile('', dummy_object=True)
            new_dataframe._segment_points.append(len(self._data))
            new_dataframe._data = pd.concat([self._data, other._data],
                                            ignore_index=True)

            return new_dataframe
        else:
            raise TypeError("One of these items isn't a Plotfile object.")

    def access(self):
        return self._data

    def get(self, key, silent=True):
        values = self._data[key]
        return values

    def hr_diagram(self, ax=None, **kwargs):
        obj = self.plot('log(T)', 'log(L)', ax=ax, **kwargs)
        xlim = obj[0].axes.get_xlim()
        if xlim[0] < xlim[1]:
            obj[0].axes.invert_xaxis()

    def kippenhahn_diagram(self, distinguish_envelopes: bool=False,
                                 legend: bool=True,
                                 x_axis: str='modelnum',
                                 ax=None,
                                 annotate: bool=True,
                                 cores_only: bool=False,
                                 **kwargs):

        err_msg = "x_axis must be collapsetime, modelnum, or age."

        assert x_axis in ['collapsetime', 'modelnum', 'age'], err_msg

        if ax is None:
            ax = plt.gca()

        TOTAL_MASS_COLOR = 'dimgrey'
        HE_CORE_MASS_COLOR = 'slateblue'
        CO_CORE_MASS_COLOR = 'crimson'
        ENVELOPE_COLOR = 'green'

        if x_axis == 'modelnum':
            X = 'timestep'
            x_label = "Model Number"
        elif x_axis == 'collapsetime':
            X = 'collapsetime'
            x_label = "Time until collapse [yr]"
        elif x_axis == 'age':
            X = 'age'
            x_label = "Age [yr]"

        if not cores_only:
            for env in range(1, 13):
                if distinguish_envelopes:
                    label = "(semi)conductive envelope"
                else:
                    label = "Conductive Envelope"

                lab = None if env < 12 else label

                transform = lambda axis,array: np.abs(array) if axis=='y' and not distinguish_envelopes else array

                self.plot(X, f'M_conv{env}', alpha=0.1,
                                             ls='',
                                             markersize=1,
                                             marker='.',
                                             label=lab,
                                             transform=transform,
                                             c=ENVELOPE_COLOR,
                                             ax=ax)

        self.plot(X, 'M', c=TOTAL_MASS_COLOR,
                          ls='-',
                          label="Total Mass",
                          ax=ax)

        self.plot(X, 'He_core',
                     c=HE_CORE_MASS_COLOR,
                     ls='-',
                     label="He core mass",
                     ax=ax) # Helium Core Mass

        self.plot(X, 'CO_core',
                     c=CO_CORE_MASS_COLOR,
                     ls='-',
                     label="CO core mass",
                     ax=ax) # CO Core Mass

        ZAMS = self.get('M').to_numpy()[0]

        if annotate:
            ax.set_title(rf"$M_{{\rm ZAMS}}={ZAMS}M_\odot$ star")
            ax.set_xlabel(x_label)
            ax.set_ylabel(r'Mass co-ordinate')

        if legend:
            ax.legend()


    def plot(self, x_axis, y_axis,
                           transform=None,
                           ax=None,
                           fix_core_masses=True,
                           **kwargs):
        """Plots the parameter given by x_axis against y_axis.

        Arguments:
            x_axis {str}           -- The x-axis to plot. Must be a key of
                                      self._data
            y_axis {str}           -- The y-axis to plot. Must be a key of
                                      self._data
            transform {callable}   -- A function to apply to each axis
            ax {Axes2D}            -- An Axes2D instance to plot on.
                                      If None, a new figure is created.
            fix_core_masses {bool} -- whether to fix the core masses such that
                                      the He core follows the envelope if the
                                      star becomes entirely stripped.
            **kwargs {variable}    -- Any keyword arguments to be passed to
                                      plt.plot.

        Notes:
            transform should be a function matching the call-signature
            transform(axis, array), where axis is 'x' or 'y' and array is the
            corresponding axis to variable `axis`. The function should return
            the transformed axis. For instance, the following function

            ```
            def do_transform(axis, array):
                if axis == 'y': return np.log10(array)
                if axis == 'x': return array + 5
            ```

            shifts the x axis by 5 and logs the y-axis content. The following
            function

            ```
            def do_transform(axis, array):
                return np.log10(array)
            ```

            logs both axes, and the following lambda:

            ```
            do_transform = lambda axis, array: array
            ```

            does nothing (and is equivalent to passing transform=None).

        Returns:
            obj {list} -- The list of Line2d objects plotted by plt.plot.
        """
        if x_axis == 'collapsetime':
            age_at_collapse = self.get('age').to_numpy()[-1]
            time_until_collapse = age_at_collapse \
                                  - self.get('age').to_numpy()

            x_arr = time_until_collapse
        else:
            x_arr = self.get(x_axis).to_numpy()

        y_arr = self.get(y_axis).to_numpy()

        if transform is not None:
            if callable(transform):
                x_arr = transform('x', x_arr)
                y_arr = transform('y', y_arr)

        if fix_core_masses and y_axis == 'He_core':
            M = self.get('M').to_numpy()
            if (y_arr != 0).any():
                mask = np.where(y_arr != 0)[0][-1]
                if mask + 1 != len(y_arr):
                    y_arr[mask:] = M[mask:]

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

    def parse_plotfile(self, fname, row, is_dummy):
        """
        Parses a plotfile.

        Arguments:
            fname    {str}  -- the filename of the plot file to be loaded
            row      {str}  -- the row that is to be loaded (obsolete here, you should specify "all", will be removed in a future version)
            is_dummy {bool} -- whether the resultant dataframe should be empty or not
        """

        c = kaitiaki.constants.PLOT_FILE_COLUMNS

        if path.exists(f'{fname}') or is_dummy:
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
                 + [12 for _ in range(18)]               # 18(1X,E12.5)
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

            if not is_dummy:
                if row == 'all':
                    df = pd.read_fwf(fname,
                                     names=c,
                                     infer_nrows=99999,
                                     colspecs='infer')
                else:
                    from file_read_backwards import FileReadBackwards

                    with FileReadBackwards(fname, encoding="utf-8") as frb:
                        line = frb.readline()

                    line = line.replace('**********', ' nan ').replace('-',' -').replace('E -','E-').split()

                    pairs = dict(zip(c, line))
                    pairs = {k: [float(v)] for k, v in pairs.items()}
                    df = pd.DataFrame.from_dict(pairs)

                status = 'loaded'
            else:
                df = pd.DataFrame(columns=c)
                status = 'dummy'
        else:
            df = None
            status = 'skipped'

        return df, status