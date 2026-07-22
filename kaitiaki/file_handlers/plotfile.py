import kaitiaki

from copy import deepcopy
from decimal import Decimal
import itertools
from os import path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import subprocess
from tqdm import tqdm
import csv

import pandas as pd


def get_last_line(self, file):
    from file_read_backwards import FileReadBackwards

    c = kaitiaki.constants.PLOT_FILE_COLUMNS

    with FileReadBackwards(file, encoding="utf-8") as frb:
        line = frb.readline()

    spec = ([6, 16]                               # I6, E16.9
            + [10 for _ in range(24)]             # 24F10.5
            + [13, 13, 13]                        # 3E13.6
            + [13 for _ in range(18)]             # 18(1X,E12.5)
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
                 engine: str = 'pandas'):

        assert row in ['all', 'last'], "row must be all or last"

        self._engine = engine

        self._data, status = self.parse_plotfile(file, row, dummy_object)
        self._segment_points = []
        self._filename = file

        if status == 'skipped':
            raise OSError("Requested file does not exist")

        self._allow_pad_age = allow_pad_age

    def basepath(self):
        return self._filename

    def __len__(self):
        return len(self._data['age'])

    def last(self):
        return self._data.iloc[-1:]

    def state(self, N):
        return self._data.iloc[N].to_numpy()

    def zams(self):
        return self._data.iloc[0:]

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
                other.pad_modelnum(self.get('timestep')[-1])

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
            values = self._data[key]
            return values.to_numpy()
        else:
            values = self._data['age']
            values = values[-1] - values

        return values.to_numpy()

    def hr_diagram(self, ax=None, limit=-1, **kwargs):
        obj = self.plot('log(T)', 'log(L)', ax=ax, limit=limit, **kwargs)

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
                           sample_every: int = 1,
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
            if distinguish_envelopes:
                label = "(semi)convective envelope"
            else:
                label = "Convective envelope"

            for env in range(1, 13):

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
                          rasterized=True,
                          sample_every=sample_every,
                          c=ENVELOPE_COLOR,
                          ax=ax)

            line = Line2D([0, 1], [0, 1], linestyle='-', color=ENVELOPE_COLOR)

        core_total = self.plot(X,
                               'M',
                               sample_every=sample_every,
                               c=TOTAL_MASS_COLOR,
                               ls='-',
                               label="Total Mass",
                               ax=ax)

        # Helium Core Mass
        core_helium = self.plot(X,
                                'He_core',
                                sample_every=sample_every,
                                c=HE_CORE_MASS_COLOR,
                                ls='-',
                                label="He core mass",
                                ax=ax)

        # CO Core Mass
        core_co = self.plot(X,
                            'CO_core',
                            sample_every=sample_every,
                            c=CO_CORE_MASS_COLOR,
                            ls='-',
                            label="CO core mass",
                            ax=ax)

        ZAMS = self.get('M')[0]

        if annotate:
            ax.set_title(rf"$M_{{\rm ZAMS}}={ZAMS}~\text{{M}}_\odot$ star")
            ax.set_xlabel(x_label)
            ax.set_ylabel(r'Mass co-ordinate / M$_\odot$')

        if legend:
            handles, labels = ax.get_legend_handles_labels()
            handles[0] = line
            labels[0] = label
            ax.legend(handles, labels, frameon=False)

    def plot(self,
             x_axis: str,
             y_axis: str,
             transform=None,
             limit=-1,
             sample_every=1,
             ax: plt.Axes = None,
             fix_core_masses: bool = True,
             **kwargs):
        """Plots the parameter given by x_axis against y_axis.

        Args:
            x_axis (str): The x-axis to plot. Must be a key of self._data
            y_axis (str): The y-axis to plot. Must be a key of self._data
            transform (callable): A function to apply to each axis
            limit (int): An integer representing how may data points to include. Unlike truncate(), this works *non-destructively*.
            sample_every (int): An integer representing how often to sample (sample_every=50 means include every 50th datapoint). Applied AFTER `limit`.
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

                obj = plt.plot(x[::sample_every], y[::sample_every], **kwargs)
                objs.append(obj)
        else:
            x = x_arr[:limit]
            y = y_arr[:limit]
            objs = [plt.plot(x[::sample_every], y[::sample_every], **kwargs)]

        return list(itertools.chain.from_iterable(objs))

    def pad_age(self, by):
        self._data['age'] += by

    def pad_modelnum(self, by):
        self._data['timestep'] += by

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
                    + [13 for _ in range(18)]            # 18(1X,E12.5)
                    + [9 for _ in range(52)])            # 52F9.5
            # The spec extends out to ~100 columns to future proof it
            # I think, so we have to truncate it here to the length of
            # what we know is in the file.
            spec = spec[:len(c)]

            if not is_dummy:
                if row == 'all':
                    import re

                    regex = r'(?<![Ee])(?<=\d)([+-]\d+)$'

                    def converter(x):
                        if x is None:
                            return np.nan

                        x = x.strip()

                        if x == '' or '*' in x:
                            return np.nan
                        return float(re.sub(regex, r'E\1', x))

                    converters = {
                        ci: converter
                        for ci in c
                    }

                    converters['timestep'] = int

                    df = pd.read_fwf(fname,
                                     names=c,
                                     widths=spec,
                                     converters=converters,
                                     # infer_nrows=99999,
                                     )
                    # else:
                    #     with open(fname, 'r') as f:
                    #         data = [[] for _ in spec]

                    #         for line in f:
                    #             pos = 0

                    #             for colnum, width in enumerate(spec):
                    #                 field = line[pos:pos + width].strip()
                    #                 if int(float(field)) == float(field):
                    #                     field = int(float(field))
                    #                 else:
                    #                     field = float(field)

                    #                 data[colnum].append(field)
                    #                 pos += width

                    #         data = {col: np.array(col_data) for col, col_data in zip(c, data)}, 'loaded'

                    #     return data, 'loaded'

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
