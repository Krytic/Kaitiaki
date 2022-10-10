from copy import deepcopy
from decimal import Decimal
from os import path

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from tqdm import tqdm

import kaitiaki

class DataFileParser:
    def __init__(self, file_pointer='data'):
        """Parses a datafile. Must be used as a context manager:

        >>> with DataFileParser('data') as dfile:
        >>>    ...

        Keyword Arguments:
            file_pointer {str} -- The location of the datafile (default: {'data'})
        """
        self._file = file_pointer

    def __enter__(self):
        class Parser():
            def __init__(self, file_pointer):
                self._datafile = file_pointer
                self._datafile_pointer = open(self._datafile, 'r+')
                self._original_contents = self._datafile_pointer.read().split("\n")
                self._contents = deepcopy(self._original_contents)

            def _check_scientific_notation(self, param):
                idx = self._get_index_of_parameter(param)
                if idx[0] == 3:
                    return 1 # everything on this line requires scientific notation to 1dp
                elif idx[0] in [17, 18]:
                    return 2 # everything on these lines requires scientific notation to 2dp
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

                """
                    God is dead. God remains dead. And we have killed him.

                    - Nietzsche
                """

                if param in lookup_table.keys():
                    return lookup_table[param]
                else:
                    raise KeyError(f"Parameter {param} not recognised.")

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

            def set(self, param, value):
                idx = self._get_index_of_parameter(param)

                value = str(value)

                num_dp_of_scientific_notation = self._check_scientific_notation(param)
                num_dp = self._determine_decimal_places(param)

                if num_dp_of_scientific_notation > 0:
                    value = f'%.{num_dp_of_scientific_notation}E' % Decimal(value)
                if num_dp > 0:
                    value = f'%.{num_dp}f' % Decimal(value)

                if idx[2] == None:
                    endpoint = len(self._contents[idx[0]])
                else:
                    endpoint = idx[2]

                length = endpoint-idx[1]
                value = value.rjust(length)

                leftpart = self._contents[idx[0]][:idx[1]]
                rightpart = self._contents[idx[0]][endpoint:]

                new_line = leftpart + value + rightpart

                self._contents[idx[0]] = new_line

                self._write_to_pointer(self._datafile_pointer, self._contents, 'replace')

            def get(self, param):
                idx = self._get_index_of_parameter(param)

                val = self._contents[idx[0]][idx[1]:idx[2]]

                num_dp = self._check_scientific_notation(param) + self._determine_decimal_places(param)

                if num_dp > 0:
                    return float(val)

                return int(val)

            def restore_backup(self):
                pass

            def set_zams_mass(self, target_mass):
                self.set('RML', target_mass)
                self.set('IML1', 9)

        self._parser = Parser(self._file)

        return self._parser

    def __exit__(self, exc_type, exc_value, traceback):
        self._parser._datafile_pointer.close()
