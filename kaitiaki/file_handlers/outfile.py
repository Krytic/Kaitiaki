import pprint

import numpy as np

import kaitiaki

def get_last_converged_model(file):
    from file_read_backwards import FileReadBackwards
    from queue import Queue

    q = []

    i = 0

    with FileReadBackwards(file, encoding="utf-8") as frb:
        for l in frb:
            while len(q) > 8:
                q.pop(0)
            q.append(l)
            if 'dt/age/MH/MHe' in l.strip():
                break

    modelblock = []

    for l in reversed(q):
        modelblock.append(l)

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
            line = line.replace('**********', ' nan ').replace('-',' -').replace('E -','E-').split()

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
            line = line.replace('**********', ' nan ').replace('-',' -').replace('E -','E-').split()
            position = line[-1]

            lines.append(line)

        self._further_information = dict()

        for i in range(len(further_info)):
            print(i)
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
                    "density, or temperature, you must specify which position "
                   f"({pos_str}) you want it at via the at keyword."))
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

class outfile:
    """Represents an *out* file object."""
    def _get_models_from_contents(self, contents, reindex_modelnums):
        models = []
        outlen = len(contents)
        model_number = 0

        for i in range(outlen):
            if 'dt/age/MH/MHe' in contents[i]:
                # we got one
                model = contents[i:i+8]
                models.append([model_number,ModelSummary(model)])

                if not reindex_modelnums:
                    models[-1][0] = int(models[-1][1].get('modelnum'))

                i += 8
                model_number += 1

        models.sort(key=lambda r: r[0])

        return models

    def __init__(self, file='out', reindex_modelnums=False):
        try:
            # this is a file path
            with open(file, 'r') as out:
                contents = out.readlines()
        except OSError:
            # this is the contents of a file
            contents = file.split("\n")

        self._models = dict(self._get_models_from_contents(contents, reindex_modelnums))

        self.__size = len(self._models)

    def get_by_modelnum(self, modelnum):
        return self._models[modelnum]

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
        return self._models[-1][1]