import kaitiaki

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
import csv


def get_last_line(self, file):
    from file_read_backwards import FileReadBackwards

    c = kaitiaki.constants.PLOT_FILE_COLUMNS

    with FileReadBackwards(file, encoding="utf-8") as frb:
        line = frb.readline()

    spec = ([6, 16]                               # I6, E16.9
            + [10 for _ in range(24)]             # 24F10.5
            + [13, 13, 13]                        # 3E13.6
            + [12 for _ in range(18)]             # 18(1X,E12.5)
            + [9 for _ in range(52)])             # 52F9.5

    for i, col in enumerate(c):
        width = spec[i]
        # TODO: Continue writing parser.


class plot:
    def __init__(self,
                 file: str = 'plot',
                 allow_pad_age: bool = True,
                 row: str = 'all',
                 dummy_object: bool = False,
                 DEBUG_load_pandas: bool = False):

        assert row in ['all', 'last'], "row must be all or last"

        if DEBUG_load_pandas is not False:
            kaitiaki.debug('warning', ('DEBUG_load_pandas is a bodge '
                                       'behaviour. It is ONLY in here to '
                                       'work as a fix for a bug that has '
                                       'been patched, though some (internal, '
                                       'unreleased) files were generated '
                                       'with it. Behaviour is untested!'))
        self._DEBUG_load_pandas = DEBUG_load_pandas

        self._data, status = self.parse_plotfile(file, row, dummy_object)
        self._segment_points = []
        self._filename = file

        if status == 'skipped':
            raise OSError("Requested file does not exist")

        self._allow_pad_age = allow_pad_age

    def __len__(self):
        return len(self._data)

    def last(self):
        return self._data.iloc[-1:]

    def zams(self):
        return self._data.iloc[0]

    def reconstruct(self, path):
        # This will break the reader...
        self._data.to_csv(path,
                          float_format='%16.3f',
                          sep=" ",
                          header=False,
                          quoting=csv.QUOTE_NONE,
                          index=False)

    def __add__(self, other):
        if isinstance(other, plot):
            if self._allow_pad_age:
                other.pad_age(self.get('age')[-1])

            # self._segment_points.append(len(self._data))
            # self._data = self._data.append(other._data)

            new_dataframe = plot('', dummy_object=True)
            new_dataframe._segment_points.append(len(self._data))
            new_dataframe._data = pd.concat([self._data, other._data],
                                            ignore_index=True)

            return new_dataframe
        else:
            raise TypeError("One of these items isn't a plot object.")

    def access(self):
        return self._data

    def get(self, key):
        if key != 'inverse_age':
            values = self._data[key].to_numpy()
            return values
        else:
            values = self._data['age'].to_numpy()
            return values[-1] - values

    def hr_diagram(self, ax=None, **kwargs):
        obj = self.plot('log(T)', 'log(L)', ax=ax, **kwargs)
        xlim = obj[0].axes.get_xlim()
        if xlim[0] < xlim[1]:
            obj[0].axes.invert_xaxis()

    def truncate(self, truncate_to):
        self._data.drop([n for n in range(truncate_to, len(self))],
                        inplace=True)

        new_dataframe = plot('', dummy_object=True)
        new_dataframe._segment_points.append(len(self._data))
        new_dataframe._data = self._data

        return new_dataframe

    def kippenhahn_diagram(self,
                           distinguish_envelopes: bool = False,
                           legend: bool = True,
                           x_axis: str = 'modelnum',
                           ax: plt.Axes = None,
                           annotate: bool = True,
                           cores_only: bool = False,
                           **kwargs):

        err_msg = "x_axis must be collapsetime, modelnum, or age."

        assert x_axis in ['collapsetime', 'modelnum', 'age'], err_msg

        if ax is None:
            ax = plt.gca()

        TOTAL_MASS_COLOR = 'dimgrey'
        HE_CORE_MASS_COLOR = 'slateblue'
        CO_CORE_MASS_COLOR = 'crimson'
        ENVELOPE_COLOR = 'green'

        match x_axis:
            case 'modelnum':
                X = 'timestep'
                x_label = "Model Number"
            case 'collapsetime':
                X = 'collapsetime'
                x_label = "Time until collapse [yr]"
            case 'age':
                X = 'age'
                x_label = "Age [yr]"

        if not cores_only:
            for env in range(1, 13):
                if distinguish_envelopes:
                    label = "(semi)convective envelope"
                else:
                    label = "Convective Envelope"

                lab = None if env < 12 else label

                if distinguish_envelopes:
                    transform = kaitiaki.utils.transforms.null
                else:
                    transform = kaitiaki.utils.transforms.absval

                self.plot(X,
                          f'M_conv{env}',
                          ls='',
                          markersize=1,
                          marker='.',
                          label=lab,
                          transform=transform,
                          c=ENVELOPE_COLOR,
                          ax=ax)

        self.plot(X,
                  'M',
                  c=TOTAL_MASS_COLOR,
                  ls='-',
                  label="Total Mass",
                  ax=ax)

        # Helium Core Mass
        self.plot(X,
                  'He_core',
                  c=HE_CORE_MASS_COLOR,
                  ls='-',
                  label="He core mass",
                  ax=ax)

        # CO Core Mass
        self.plot(X,
                  'CO_core',
                  c=CO_CORE_MASS_COLOR,
                  ls='-',
                  label="CO core mass",
                  ax=ax)

        ZAMS = self.get('M')[0]

        if annotate:
            ax.set_title(rf"$M_{{\rm ZAMS}}={ZAMS}M_\odot$ star")
            ax.set_xlabel(x_label)
            ax.set_ylabel(r'Mass co-ordinate')

        if legend:
            ax.legend()

    def plot(self,
             x_axis: str,
             y_axis: str,
             transform=None,
             ax: plt.Axes = None,
             fix_core_masses: bool = True,
             **kwargs):
        """Plots the parameter given by x_axis against y_axis.

        Args:
            x_axis (str): The x-axis to plot. Must be a key of self._data
            y_axis (str): The y-axis to plot. Must be a key of self._data
            transform (callable): A function to apply to each axis
            ax (Axes2D): An Axes2D instance to plot on. If None, a new figure is created.
            fix_core_masses (bool): Whether to fix the core masses such that the He core follows the envelope if the star becomes entirely stripped.
            **kwargs (variable): Any keyword arguments to be passed to plt.plot.

        Returns:
            list: The list of Line2d objects plotted by plt.plot.

        Notes:
            transform should be a function matching the call-signature :code:`transform(axis, array)`, where axis is 'x' or 'y' and array is the corresponding axis to variable :code:`axis`. The function should return the transformed axis. For instance, the following function::

                def do_transform(axis, array):
                    if axis == 'y': return np.log10(array)
                    if axis == 'x': return array + 5

            shifts the x axis by 5 and logs the y-axis content. The following
            function::

                def do_transform(axis, array):
                    return np.log10(array)

            logs both axes, and the following lambda:::

                do_transform = lambda axis, array: array

            does nothing (and is equivalent to passing transform=None).
        """

        if x_axis == 'collapsetime':
            age_at_collapse = self.get('age')[-1]
            current_age = self.get('age')

            time_until_collapse = age_at_collapse - current_age

            x_arr = time_until_collapse
        else:
            x_arr = self.get(x_axis)

        y_arr = self.get(y_axis)

        if transform is not None:
            if callable(transform):
                x_arr = transform('x', x_arr)
                y_arr = transform('y', y_arr)

        if fix_core_masses and y_axis == 'He_core':
            M = self.get('M')
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

    def parse_plotfile(self, fname: str = 'plot',
                       row: str = 'all',
                       is_dummy: bool = False):
        """Parses a plotfile.

        Args:
            fname (str): The plot file to load
            row (str): The row to load (obsolute -- always pass "all")
            is_dummy (bool): Whether the dataframe should be empty or not

        Returns:
            [type]: [description]
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
            spec = ([6, 16]                               # I6, E16.9
                    + [10 for _ in range(24)]            # 24F10.5
                    + [13, 13, 13]                       # 3E13.6
                    + [12 for _ in range(18)]            # 18(1X,E12.5)
                    + [9 for _ in range(52)])            # 52F9.5
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
                    if self._DEBUG_load_pandas:
                        # Don't use! This is to work around a patched bug
                        df = pd.read_csv(fname,
                                         names=['index', *c],
                                         header=0,
                                         sep=r'\s+',
                                         index_col=0,
                                         dtype=float)
                    else:
                        df = pd.read_fwf(fname,
                                         names=c,
                                         # widths=spec
                                         infer_nrows=99999,
                                         )
                else:
                    from file_read_backwards import FileReadBackwards

                    with FileReadBackwards(fname, encoding="utf-8") as frb:
                        line = frb.readline()

                    line = (line.replace('**********', ' nan ')
                                .replace('-', ' -')
                                .replace('E -', 'E-')
                                .split())

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


def plot2(file: str = 'plot2', *args, **kwargs):
    return plot(file, *args, **kwargs)
