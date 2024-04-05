import pprint
import tempfile

import numpy as np

import kaitiaki

from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class DeprecationWarning(Exception):
    pass


def get_last_converged_model(file):
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

    return ModelSummary(modelblock)


class ModelSummary:
    def _parse_block(self, model):
        self._datapoints = dict()
        self._originalformat = model
        model_num = float(model[0].split()[0])

        lines = []

        column_headers = model[0].split()[1:]
        self._column_headers = column_headers

        for lineno in range(1, 8):
            line = model[lineno]
            line = (line.replace('**********', ' nan ')
                        .replace('-', ' -')
                        .replace('********', ' nan ')
                        .replace('E -', 'E-')
                        .split())

            lines.append(line)

            # Now we extract abundances, degeneracy, density, and temperature
            # at the three points of interest: center, surface, and point
            # of maximum temperature
            # This is only valid for the first three lines:
            if lineno < 4:
                position = line[-1]
                data = dict()
                for colno in range(6, 16):
                    data[column_headers[colno-1]] = float(line[colno])

                self._datapoints[position] = data

        self._information = dict()

        informations = [['mass', 'dt', 'tn', 'P', 'LH', 'Lth'],
                        ['H_env', 'age', 'tKH', 'rlf', 'LHe', 'Lnu'],
                        ['He_env', 'MH', 'MHe', 'Mb', 'dM', 'LC', 'm'],
                        ['MC', 'conrad1', 'conrad2', 'conrad3', 'conrad4',
                         'conrad5', 'conrad6', 'conrad7', 'conrad8', 'conrad9',
                         'conrad10', 'conrad11', 'conrad12', 'logmr2k2', 'k2',
                         'logR', 'logL', 'nmod2'],
                        ['k_H_bound', 'k_He_bound', 'k_conrad1', 'k_conrad2',
                         'k_conrad3', 'k_conrad4', 'k_conrad5', 'k_conrad6',
                         'k_conrad7', 'k_conrad8', 'k_conrad9', 'k_conrad10',
                         'k_conrad11', 'k_conrad12'],
                        ['M_maxH', 'M_maxHe', 'M_maxC', 'M_burn1', 'M_burn2',
                         'M_burn3', 'M_burn4', 'M_burn5', 'M_burn6', 'M_burn7',
                         'M_burn8', 'M_burn9', 'M_burn10', 'M_burn11',
                         'M_burn12'],
                        ['k_maxH', 'k_maxHe', 'k_maxC', 'k_burn1', 'k_burn2',
                         'k_burn3', 'k_burn4', 'k_burn5', 'k_burn6', 'k_burn7',
                         'k_burn8', 'k_burn9', 'k_burn10', 'k_burn11',
                         'k_burn12']]

        for i in range(7):
            for j in range(len(informations[i])):
                label = informations[i][j]

                self._information[label] = float(lines[i][j])

        lines = []

        self._further_information = dict()

        return

        further_info = []

        for lineno in range(4, 4+len(further_info)):
            line = model[lineno]
            line = (line.replace('**********', ' nan ')
                        .replace('-', ' -')
                        .replace('E -', 'E-')
                        .split())

            position = line[-1]

            lines.append(line)

        self._further_information = dict()

        for i in range(len(further_info)):
            for j in range(len(further_info[i])):
                print(f"    {j}")
                label = further_info[i][j]

                self._further_information[label] = float(lines[i][j])

    def __init__(self, model):
        i = 0
        modelnum = model[i].strip().split()[0]
        MHe = model[i+2].strip().split()[0]
        mass = model[i+1].strip().split()[0]

        self._modelnum = float(modelnum)
        self._mass = float(mass)

        self._parse_block(model)

    def nan_num(self):
        attrs = self.all_attributes()

        nans = 0

        for k, v in attrs.items():
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    # max depth is 2
                    if np.isnan(v2):
                        nans += 1
            else:
                if np.isnan(v):
                    nans += 1

        return nans

    def original_format(self):
        return "\n".join(self._originalformat)

    def get(self, param, at=None):
        if param == 'modelnum':
            return self._modelnum
        elif param in self._column_headers:
            if at not in self._datapoints.keys():
                pos_str = "/".join(self._datapoints.keys())
                raise ValueError(("If you request an abundance, degeneracy, "
                                  "density, or temperature, you must specify "
                                 f"which position ({pos_str}) you want it at "
                                  "via the `at` keyword."))
            else:
                return self._datapoints[at][param]
        elif param in self._information.keys():
            return self._information[param]
        elif param in self._further_information.keys():
            return self._further_information[param]
        else:
            raise KeyError(f'{param} not found')

    def all_attributes(self):
        return {**self._datapoints, **self._information}

    def __str__(self):
        attrs = self.all_attributes()

        pp = pprint.PrettyPrinter(indent=4)
        return pp.pformat(attrs)


class ModelProfile:
    def __init__(self, array):
        self._array = array


class outfile:
    """Represents an *out* file object."""
    def _get_models_from_contents(self, contents, reindex_modelnums):
        models = []
        outlen = len(contents)
        model_number = 0

        for i in range(outlen):
            if 'dt/age/MH/MHe' in contents[i]:
                # we got a model summary
                model = contents[i:i+8]
                models.append([model_number, ModelSummary(model)])

                if not reindex_modelnums:
                    models[-1][0] = int(models[-1][1].get('modelnum'))

                i += 8
                model_number += 1

        models.sort(key=lambda r: r[0])

        return models

    def _get_profiles_from_contents(self, contents, K=199):
        profiles = []

        for line in contents:
            if len(line.split()) < 2:
                continue
            elif line.split()[0] == 'K':
                dtypes = [('k', 'int')] + \
                         [(name, 'float') for name in line.split()[1:]]
                data = [contents[i] for i in range(K+(K-1)//10)]
                data = [tuple(map(float, row.split()))
                        for row in data if row != '\n']
                profiles.append(np.array(data, dtype=dtypes))

        return profiles

    def _tomso(self, filename):
        """Reads a STARS `out` file and returns (part of) the summaries and
        profiles in two structured arrays.

        Parameters
        ----------
        filename: str
            Filename of the STARS output file to load.

        Returns
        -------
        summaries: 2-d structured array
            Summaries of each model in the run, similar to MESA's
            `history` files.

        profiles: 3-d structured array
            Model profiles produced at regular intervals during the run.
            The first index of the array is the profile number.
        """
        get_data = True
        if filename.endswith('out2'):
            with open(filename[:-1], 'r') as f:
                line = f.readline()

                summaries = []
                profiles = []

                data = [line.rstrip('\n')]

                while (line := f.readline()).strip() != '':
                    # read in data
                    data.append(line.rstrip('\n'))

                _fp = tempfile.NamedTemporaryFile('w+')

                _fp.write("\n".join(data))
                _fp.flush()
                _fp.seek(0)

                self.data_block = kaitiaki.file.data(_fp.name)
                self.__fp = _fp

                K = self.data_block.get('NM2')

                get_data = False

        with open(filename, 'r') as f:
            if get_data:
                line = f.readline()

                summaries = []
                profiles = []

                data = [line.rstrip('\n')]

                while (line := f.readline()).strip() != '':
                    # read in data
                    data.append(line.rstrip('\n'))

                _fp = tempfile.NamedTemporaryFile('w+')

                _fp.write("\n".join(data))
                _fp.flush()
                _fp.seek(0)

                self.data_block = kaitiaki.file.data(_fp.name)
                self.__fp = _fp

                K = self.data_block.get('NM2')

            headers = []

            while (line := f.readline()):
                if len(line.split()) < 2:
                    continue
                elif line.split()[0] == 'K':
                    # found a profile
                    dtypes = [('k', 'int')] + \
                             [(name, 'float') for name in line.split()[1:]]
                    data = [f.readline() for i in range(K+(K-1)//10)]
                    data = [tuple(map(float, row.split()))
                            for row in data if row != '\n']
                    profiles.append(np.array(data, dtype=dtypes))

                    if dtypes not in headers:
                        headers.append(dtypes)

        # Handle "pagination".
        # If NWRT3 = 2, for example, that means that every profile comes in
        # a set of 2, i.e., profiles[0] and profiles[1] represent the same
        # timestep, just different data "pages".
        # So we need to create a list of tuples:
        # (profiles[0], profiles[1], ... profiles[NWRT3-1], profiles[NWRT3].
        def partition(L, N): return list(zip(*([iter(L)] * N)))
        # Function above adapted from thefourtheye on StackOverflow:
        # https://stackoverflow.com/questions/23286254/how-to-convert-a-list-to-a-list-of-tuples

        NWRT3 = self.data_block.get('NWRT3')

        self.pages_per_profile = NWRT3
        self.profile_headers = np.vstack((headers,))
        self.n_profiles = len(profiles) / NWRT3

        profiles = partition(profiles, NWRT3)

        return np.vstack((profiles,))

    def __init__(self, file='out', reindex_modelnums=False):
        try:
            # this is a file path
            with open(file, 'r') as out:
                contents = out.readlines()
        except OSError:
            # this is the contents of a file
            raise DeprecationWarning(('Passing the file contents as a '
                                      'parameter to kaitiaki.file.out is '
                                      'deprecated. Please pass a path to '
                                      'an outfile.'))

        models = self._get_models_from_contents(contents, reindex_modelnums)

        self._models = dict(models)

        self._profiles = self._tomso(file)

        self.__size = len(self._models)

        self.__detailed_idx = 0

    def __del__(self):
        self.__fp.close()

    def get_by_modelnum(self, modelnum):
        return self._models[modelnum]

    def get(self, param, at=None):
        ret = np.array([])

        for model in self:
            ret = np.append(ret, model.get(param, at))

        return ret

    def detailed_profile(self, idx):
        return self._profiles[idx]

    def detailed_profile_history(self):
        raise NotImplementedError()
        self.__detailed_idx += 1
        yield self.detailed_profile(self.__detailed_idx-1)

    def __iter__(self):
        self.__counter = 0
        return self

    def __next__(self):
        if self.__counter < self.__size-1:
            self.__counter += 1
            return self.get_by_modelnum(self.__counter-1)
        else:
            raise StopIteration

    def __len__(self):
        return self.__size

    def last(self):
        key = max(self._models.keys())
        return self._models[key]

    def find_column(self, param: str):
        NWRT3 = self.pages_per_profile
        headers = self.profile_headers

        column_location = None
        column_pagenum = None

        for pagenum in range(NWRT3):
            if param in headers[pagenum][:, 0]:
                cond = headers[pagenum][:, 0] == param
                column_location = np.where(cond)[0][0]
                column_pagenum = pagenum
                break
        else:
            raise KeyError(f"Parameter {param} not found in the out file.")

        return column_location, column_pagenum

    def build_network(self, param: str):
        N = int(self.n_profiles)

        NM2 = self.data_block.get('NM2')
        NWRT1 = self.data_block.get('NWRT1')

        meshpoints = np.arange(1, NM2+1)

        results = np.zeros((N, NM2))

        network_index = np.arange(N) * NWRT1

        column_location, column_pagenum = self.find_column(param)

        if column_location is not None:
            for idx in range(N):
                profile = self.detailed_profile(idx)

                page = profile[column_pagenum]

                column_data = [row[column_location] for row in page]

                model_number = idx * NWRT1

                results[idx, :] = column_data

        return network_index, meshpoints, results

    def visualise_parameter_network(self,
                                    param: str,
                                    mode: str = '2d',
                                    nrows: int = 1,
                                    ncols: int = 1,
                                    fidelity: int = 100,
                                    span: tuple = (None, None),
                                    make_cbar: bool = True,
                                    figsize=None):

        assert mode.lower() in ['2d', '3d'], 'mode must be 2d or 3d.'
        assert nrows > 0, 'nrows must be positive.'
        assert ncols > 0, 'ncols must be positive.'

        N = int(self.n_profiles)
        NWRT1 = self.data_block.get('NWRT1')

        mode = mode.lower()

        network_index, meshpoints, results = self.build_network(param)

        X = network_index
        Y = meshpoints
        XX, YY = np.meshgrid(X, Y, indexing='ij')

        fig = plt.figure(figsize=figsize)

        if mode == '2d':
            ax = fig.add_subplot(nrows, ncols, 1)

            if span == (None, None):
                span = (results.min(), results.max())

            levels = np.linspace(span[0], span[1], fidelity)

            cmap = ax.contourf(XX, YY, results, cmap=cm.bwr,
                               levels=levels)

            ax.set_xlabel('Model Number')
            ax.set_ylabel('Meshpoint Number')

            if make_cbar:
                cbar = fig.colorbar(cmap)
                cbar.set_label(param)
        else:
            ax = fig.add_subplot(nrows, ncols, 1, projection='3d')

            surf = ax.plot_surface(XX, YY, results, cmap=cm.bwr,
                                   linewidth=0, antialiased=False)

            ax.set_xlabel('Model Number')
            ax.set_ylabel('Meshpoint Number')
            ax.set_zlabel(param)

        return fig


class outfile2(outfile):
    def __init__(self, file_location="out2", reindex_modelnums=False):
        super().__init__(file_location, reindex_modelnums)
