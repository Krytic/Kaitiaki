from copy import deepcopy
from decimal import Decimal
from os import path
import pickle
import tempfile

import colorama
from colorama import Fore, Back, Style
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tabulate import tabulate
from tqdm import tqdm

import kaitiaki


class DataFileParser:
    def __init__(self, file_pointer='data'):
        """Parses a datafile. Must be used as a context manager:

        >>> with DataFileParser('data') as dfile:
        >>>    ...

        If not used as a context manager, only the get(), as_dict(), and
        comparison dunder methods are available.

        Keyword Arguments:
            file_pointer {str} -- The location of the datafile
                                  (default: {'data'})
        """
        self._file = file_pointer

    def get(self, param):
        with self as data:
            return data.get(param)

    def elucidate(self, *args, **kwargs):
        with self as data:
            return data.elucidate(*args, **kwargs)

    def to_pickle(self, path: str = "configuration.params"):
        params = self.as_dict()

        with open(path, 'wb') as handle:
            pickle.dump(params, handle)

    def explain(self, param):
        val = self.get(param)
        if param in kaitiaki.constants.disambiguable:
            meaning = kaitiaki.constants.disambiguable[param]

            lookup = meaning['options'].items()
            valid = [f"    {key}: {value}" for key, value in lookup]
            valid = '\n'.join(valid)

            disambiguation = f"""
===================
== Disambiguator ==
===================

Selected: {param.upper()} ({meaning['description']})
Valid Values:
{valid}

Current Value: {val} ({meaning['options'][val]})
            """

            print(disambiguation)

    def show_in_file(self, param):
        param = param.lower()
        if param not in kaitiaki.constants.dfile_struct.keys():
            raise KeyError('Invalid parameter')

        with self as data:
            # Now we're cooking
            colorama.init()
            row_loc, start, finish = kaitiaki.constants.dfile_struct[param]
            substr = data._contents[row_loc][start:finish]

            infos = kaitiaki.constants.dfile_struct.values()
            j = max([info[0] for info in infos])  # Get the max row number

            def hl(string):
                return Fore.GREEN + Style.BRIGHT + string + Style.RESET_ALL

            for i, row in enumerate(data._contents):
                if i > j:
                    break
                if i != row_loc:
                    print(row)
                else:
                    if finish is None:
                        end = ''
                    else:
                        end = row[finish:]
                    print(row[:start] + hl(substr) + end)

    def as_dict(self):
        keys = list(kaitiaki.constants.dfile_struct.keys())

        return self.get(keys)

    def __eq__(self, other):
        my_values = self.as_dict()
        other_val = other.as_dict()

        return my_values == other_val

    def changes_from_base(self):
        base_dfile_contents = kaitiaki.load_file('data.bak')

        # I hate this implementation -- better to get a path to data.bak
        # and then load that here...
        _fp = tempfile.NamedTemporaryFile('w+')

        _fp.write("\n".join(base_dfile_contents))
        _fp.flush()
        _fp.seek(0)

        print(base_dfile_contents)

        base_dfile_object = kaitiaki.file.data(_fp.name)

        mismatches = self.compare(base_dfile_object)

        for key, val in mismatches.items():
            mismatches[key] = {
                'current': val['self'],
                'base': val['other']
            }

        return mismatches

    def __str__(self):
        with self as data:
            return '\n'.join(data._contents)

    def compare(self, other, tag_names=('other', 'self')):
        mismatches = dict()

        for key in kaitiaki.constants.dfile_struct.keys():
            other_val = other.get(key)
            my_val = self.get(key)

            if other_val != my_val:
                mismatches[key] = {tag_names[0]: other_val,
                                   tag_names[1]: my_val}

        return mismatches

    def __enter__(self):
        class Parser():
            def __init__(self, file_pointer):
                self._datafile = file_pointer
                self._datafile_pointer = open(self._datafile, 'r+')
                self._original_contents = (self._datafile_pointer.read()
                                                                 .split("\n"))
                self._contents = deepcopy(self._original_contents)

            def __str__(self):
                return "\n".join(self._contents)

            def set_from_pickle(self, path):
                with open(path, 'rb') as handle:
                    params = pickle.load(handle)

                self.set(params)

            def _check_scientific_notation(self, param):
                idx = self._get_index_of_parameter(param)
                if idx[0] == 3:
                    # This line requires scientific notation to 1dp
                    return 1
                elif idx[0] in [17, 18]:
                    # These lines require scientific notation to 2dp
                    return 2
                else:
                    if param.lower() in ['ct8', 'ct9', 'ct10']:
                        return 1
                    if param.lower() == 'zs':
                        return 2
                    if param.lower() in ['vrot1', 'vrot2']:
                        return 2
                    if param.lower() in ['facsgmin', 'sgthfac']:
                        return 2
                    if param.lower() in ['hkh', 'gff']:
                        return 2

                return 0

            def _determine_decimal_places(self, param):
                idx = self._get_index_of_parameter(param)
                param = param.lower()
                if idx[0] == 16 and param != 'zs':
                    return 3
                if param in ['fmac', 'fam']:
                    return 2
                if param in ['trc1', 'trc2']:
                    return 1
                if param == 'mwts':
                    return 2
                if idx[0] == 15 and param not in ['ct10', 'ct9', 'ct8']:
                    return 2
                if param == 'alphace':
                    return 1

                return 0

            def _get_index_of_parameter(self, param):
                """
                                                  .

                                                   .
                                         /^\     .
                                    /\   "V"
                                   /__\   I      O  o
                                  //..\\  I     .S
                                  \].`[/  I
                                  /l\/j\  (]    .  O
                                 /. ~~ ,\/I          .
                                 \\L__j^\/I       o
                                  \/--v}  I     o   .
                                  |    |  I   _________
                                  |    |  I c(`       ')o
                                  |    l  I   \.     ,/
                                _/j  L l\_!  _//^---^\\_
                                     Here be wizard.
                """
                param = param.lower()
                # Structure of the lookup table:
                # 3-tuple: (line number, starting point, ending point)
                lookup_table = kaitiaki.constants.dfile_struct

                if param in lookup_table.keys():
                    return lookup_table[param]
                else:
                    raise KeyError(f"Parameter {param} not recognised.")

            def _format_elucidate(self, datadict, fmt, splice_into=2):
                if fmt == 'plaintext':
                    res = ''
                    for key, val in datadict.items():
                        res += f"[{key}]: {val}\n"
                else:
                    # Gotta be a cleaner way to do this...
                    def format(item):
                        item = str(item)
                        item = item.lower()

                        if 'e' in item:
                            base, mantissa = item.split('e')

                            def num_is(a, b):
                                if float(a) == int(float(a)):
                                    if float(a) == float(b):
                                        return True
                                return False

                            if num_is(base, 0):
                                return 0
                            if num_is(base, 1):
                                prefix = ''
                            else:
                                prefix = f"{base}\\times"
                            mantissa = int(mantissa)
                            item = f'${prefix}10^{{{mantissa}}}$'
                        else:
                            if float(item) > 1e3 or float(item) < 1e-3:
                                item = format(f"{Decimal(item):.3E}")

                        return item

                    table = [[k, format(v), ''] for k, v in datadict.items()]

                    new_shape = [int(len(table)//splice_into), 3*splice_into]

                    if new_shape[0]*new_shape[1] < len(table)*3:
                        new_shape[0] += 1

                        print(f"{len(table)*3=}, target {new_shape[0]*new_shape[1]}")

                        print(new_shape)

                        while len(table)*3 < new_shape[0]*new_shape[1]:
                            print(f"{len(table)=}")
                            table.append(['-', '-', '-'])

                        print(table)

                    table = np.reshape(np.array(table), new_shape)
                    res = tabulate(table,
                                   headers=['Parameter', 'Value', 'Comment']*splice_into,
                                   tablefmt='latex_raw',
                                   colalign=["right", "left", "left"]*splice_into)
                return res

            def elucidate(self, file, elucidate_only=None,
                                      fmt='plaintext',
                                      splice_into=2):
                error = "`fmt` must be tex or plaintext"
                assert fmt in ['plaintext', 'tex'], error

                with open(file, 'w') as f:
                    if elucidate_only is None:
                        elucidate_only = kaitiaki.constants.dfile_struct.keys()

                    lines = dict()

                    for key in elucidate_only:
                        # ISX16-ISX18 fail for... some reason.
                        if key.lower().startswith('isx'): continue
                        if key.lower().startswith('nuc'): continue
                        if key.lower().startswith('evo'): continue

                        lines[key.upper()] = self.get(key)

                    elucidation = self._format_elucidate(lines,
                                                         fmt,
                                                         splice_into)

                    f.writelines(elucidation)

            def _write_to_pointer(self, pointer, data, mode):
                if mode == 'replace':
                    pointer.seek(0)
                pointer.write("\n".join(data))

            def make_backup(self):
                with open(self._datafile + ".bak", 'w') as file:
                    self._write_to_pointer(file, self._original_contents, 'w')

            def backup_if_not_exists(self):
                if not path.exists(self._datafile + ".bak"):
                    self.make_backup()

            def _setitem(self, param, value):
                idx = self._get_index_of_parameter(param)

                value = str(value)

                if param.lower() == 'zs':
                    if value[0] == 'z':
                        if value[1:3] != 'em':
                            value = '0.' + value[1:]
                        else:
                            value = str(10**(-int(value[3:])))

                scientific_notation_dp = self._check_scientific_notation(param)
                num_dp = self._determine_decimal_places(param)

                if scientific_notation_dp > 0:
                    value = f'%.{scientific_notation_dp}E' % Decimal(value)
                if num_dp > 0:
                    value = f'%.{num_dp}f' % Decimal(value)

                if idx[2] is None:
                    endpoint = len(self._contents[idx[0]])
                else:
                    endpoint = idx[2]

                length = endpoint-idx[1]
                value = value.rjust(length)

                leftpart = self._contents[idx[0]][:idx[1]]
                rightpart = self._contents[idx[0]][endpoint:]

                new_line = leftpart + value + rightpart

                self._contents[idx[0]] = new_line

                self._write_to_pointer(self._datafile_pointer,
                                       self._contents,
                                       'replace')

            def set(self, *args):
                if len(args) == 1:
                    try:
                        for key, value in args[0].items():
                            self._setitem(key, value)
                    except AttributeError:
                        raise TypeError(("If only one argument is passed to "
                                         "set(), it should be a dictionary "
                                         "of key-value pairs."))
                elif len(args) == 2:
                    self._setitem(*args)
                else:
                    raise TypeError(("Malformed arguments to set(). "
                                     "I expected a dictionary of key-value "
                                     "pairs, or a key and a value."))

            def _getitem(self, param):
                idx = self._get_index_of_parameter(param)

                val = self._contents[idx[0]][idx[1]:idx[2]]

                decimal_places = self._determine_decimal_places(param)
                scientific_notation = self._check_scientific_notation(param)

                num_dp = scientific_notation + decimal_places

                if num_dp > 0:
                    return float(val)

                try:
                    return int(val)
                except ValueError as e:
                    raise ValueError(f"Failed to parse parameter {param=}")

            def get(self, param):
                if isinstance(param, list):
                    ret_array = dict()

                    for item in param:
                        ret_array[item] = (self._getitem(item))

                    return ret_array

                return self._getitem(param)

            def restore_backup(self):
                pass

            def set_zams_mass(self, target_mass):
                self.set('RML', target_mass)
                self.set('IML1', 9)

        self._parser = Parser(self._file)

        return self._parser

    def __exit__(self, exc_type, exc_value, traceback):
        self._parser._datafile_pointer.close()
