import os
from io import StringIO

import kaitiaki
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Model:
    def __init__(self, model_df, delta_df, metadata, raw):
        self.__data = {
            'model': model_df,
            'delta': delta_df,
            'metadata': metadata,
            'raw': raw
            }

    # Backwards compatibility
    def __getitem__(self, item):
        if item in self.__data.keys():
            return self.__data[item]

    def text(self):
        return self.__data['raw']


class ModelFile:
    def __init__(self, file_location):
        if not os.path.exists(file_location):
            raise FileNotFoundError(f'{file_location} does not exist.')

        self.__models = []

        with open(file_location, 'r') as f:
            contents = f.readlines()

        first_line = (contents[0].replace('-', ' -')
                                 .replace('E -', 'E-')
                                 .split())

        nmesh = int(first_line[6])

        lines_per_model = 2 * nmesh + 1

        self.__NM2 = nmesh

        self.__length = 0

        for i in range(0, len(contents), lines_per_model):
            data = contents[i:i+lines_per_model]

            header = data[0].strip().split()

            model = ("".join(data[1:nmesh+1])
                       .replace('-', ' -')
                       .replace('E -', 'E-'))

            deltas = ("".join(data[nmesh+1:])
                        .replace('-', ' -')
                        .replace('E -', 'E-'))

            model_df = pd.read_csv(StringIO(model),
                                   names=['log(f)', 'log(T)', 'X_O',
                                          'log(m)', 'X_H', 'dQ/dK',
                                          'log(r)', 'L', 'X_He',
                                          'X_C', 'X_Ne', 'X_N',
                                          'H_orb', 'H_spin', 'X_3He'],
                                   sep=r"\s+")

            delta_df = pd.read_csv(StringIO(model),
                                   names=['log(f)', 'log(T)', 'X_O',
                                          'log(m)', 'X_H', 'dQ/dK',
                                          'log(r)', 'L', 'X_He',
                                          'X_C', 'X_Ne', 'X_N',
                                          'H_orb', 'H_spin', 'X_3He'],
                                   sep=r"\s+")

            first_line = (contents[i].replace('-', ' -')
                                     .replace('E -', 'E-')
                                     .split())

            keys = ['mass', 'dt', 'age', 'period', 'total_mass',
                    'energy_generation', 'NM2', 'N_models_target',
                    'NM_start', 'ISTAR', 'H_pressure', 'He_pressure']

            model_df.index += 1
            delta_df.index += 1

            if len(first_line) != len(keys):
                # For some reason, there's sometimes no pressure terms.
                keys = keys[:-2]

            metadata = dict(zip(keys, first_line))

            model = Model(model_df, delta_df, metadata, "".join(data))

            self.__models.append(model)

            self.__length += 1

    def __len__(self):
        return self.__length

    def request_indexes(self):
        idxs = np.array([])
        for i in range(len(self)):
            idxs = np.append(idxs, self.get(i)['metadata']['NM_start'])

        return idxs.astype(int)

    def get(self, index):
        return self.__models[index]

    def get_by_modelnum(self, modelnum):
        idxs = self.request_indexes()
        nbr = np.argwhere(idxs == modelnum)[0][0]
        return self.get(nbr)

    def plot(self,
             parameter,
             timestep_index=0,
             transformation=None,
             labelname=""):

        if transformation is None:
            transformation = kaitiaki.utils.transforms.null

        model = self.get(timestep_index)
        df = model['model']
        NM2 = int(model['metadata']['NM2'])

        values = transformation('y', df[parameter].to_numpy()[::-1])
        # Reversing index because the dataframes go 1 -> 199
        # (i.e., surface to center), but we want to go center -> surface.

        kaitiaki.augments.radial_plot(values, NM2, labelname)

        plt.show()


class modin(ModelFile):
    def __init__(self, file_location='modin'):
        super().__init__(file_location)


class modout(ModelFile):
    def __init__(self, file_location='modout'):
        super().__init__(file_location)


class modin2(ModelFile):
    def __init__(self, file_location='modin2'):
        super().__init__(file_location)


class modout2(ModelFile):
    def __init__(self, file_location='modout2'):
        super().__init__(file_location)
