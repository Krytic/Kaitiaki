import pprint

import numpy as np

import kaitiaki

def get_last_converged_model(self, file):
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

    return kaitiaki.out.ModelSummary(modelblock)

class ModelSummary:
    def _parse_block(self, model):
        self._datapoints = dict()

        model_num = float(model[0].split()[0])

        lines = []

        column_headers = model[0].split()[1:]
        self._column_headers = column_headers

        for lineno in range(1, 4):
            line = model[lineno]
            line = line.replace('**********', ' nan ').replace('-',' -').replace('E -','E-').split()
            position = line[-1]

            lines.append(line)

            # Now we extract abundances, degeneracy, density, and temperature
            # at the three points of interest: center, surface, and point
            # of maximum temperature
            data = dict()

            for colno in range(6, 16):
                data[column_headers[colno-1]] = float(line[colno])

            self._datapoints[position] = data

        self._information = dict()

        informations = [['mass', 'dt', 'tn', 'P', 'LH', 'Lth'],
                        ['H_env', 'age', 'tKH', 'rlf', 'LHe', 'Lnu'],
                        ['He_env', 'MH', 'MHe', 'Mb', 'dM', 'LC', 'm']]

        for i in range(3):
            for j in range(len(informations[i])):
                label = informations[i][j]

                self._information[label] = float(lines[i][j])

    def __init__(self, model):
        i = 0
        modelnum = model[i].strip().split()[0]
        MHe = model[i+2].strip().split()[0]
        mass = model[i+1].strip().split()[0]

        self._modelnum = float(modelnum)
        self._mass = float(mass)

        self._model_data = self._parse_block(model)

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
        else:
            return self._information[param]

    def all_attributes(self):
        return {**self._datapoints, **self._information}

    def __str__(self):
        attrs = self.all_attributes()

        pp = pprint.PrettyPrinter(indent=4)
        return pp.pformat(attrs)

class outfile:
    """Represents an *out* file object."""
    def _get_models_from_contents(self, contents):
        models = []
        outlen = len(contents)
        for i in range(outlen):
            if 'dt/age/MH/MHe' in contents[i]:
                # we got one
                model = contents[i:i+8]
                models.append([0,ModelSummary(model)])
                models[-1][0] = models[-1][1].get('modelnum')
                i += 8

        models.sort(key=lambda r: r[0])

        return models

    def __init__(self, file='out'):
        try:
            # this is a file path
            with open(file, 'r') as out:
                contents = out.readlines()
        except OSError:
            # this is the contents of a file
            contents = file.split("\n")

        self._models = self._get_models_from_contents(contents)

    def get_by_modelnum(self, modelnum):
        return dict(self._models)[modelnum]

    def last(self):
        return self._models[-1][1]