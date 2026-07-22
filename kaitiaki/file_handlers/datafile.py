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
    """Representation of a datafile object

    IMPORTANT NOTE:
        The parser only allows read-only operations if not used as a
        context manager. This is intentional design to limit the risk
        of conflicting write operations, and to ensure that the file is
        correctly closed after it is written to.
    """
    def __init__(self, file_pointer='data'):
        """Parses a datafile. Must be used as a context manager:

        >>> with DataFileParser('data') as dfile:
        >>>    ...

        If not used as a context manager, only the get(), as_dict(), and
        comparison dunder methods are available.

        Keyword Arguments:
            file_pointer (str): The location of the datafile
                                (default: 'data')
        """
        self._file = file_pointer

    def get(self, param: str):
        """Retrieve a parameter from the datafile.

        Retrieves a parameter by name from the datafile, and returns the value.
        Syntactic sugar for the Parser::data subclass.

        Args:
            param (str): The parameter we wish to retrieve. Must be a named
                         parameter that exists in the datafile.

        Returns:
            mixed: The parameter requested. This is returned as a float if the
                   parameter is float-like in the datafile, and an int otherwise.

        Raises:
            ValueError: if the parameter does not exist in the datafile.
        """
        with self as data:
            return data.get(param)

    def elucidate(self, *args, **kwargs):
        """Generate a report of all of the parameters in the data file as
        key-value pairs.

        Creates a file with key: value pairs of parameters. This file is
        intended to be human-readable rather than machine-readable (you
        should be using the raw datafile for machine-readable methods
        anyways).

        Args:
            file (str): The file to output the elucidation to.

        Keyword Args:
            eludicate_only (list): A list of parameters to elucidate. If
                                   None, then all parameters in the datafile
                                   will be elucidated. (default: None)

            fmt (str): The file format to elucidate (either `tex` or
                       `plaintext`). `tex` creates a booktabs table to be
                       pasted into a TeX document, whereas `plaintext` is
                       just a list of key-value pairs.

            splice_into (int): An integer representing how many column-pairs
                               to splice the list into if we are printing a
                               tex'd table. For instance, splice_int0=2 gives
                               a table like:
                               Parameter | Value | Comment | Parameter | Value | Comment |
        """
        with self as data:
            return data.elucidate(*args, **kwargs)

    def to_pickle(self, path: str = "configuration.params"):
        """Writes the datafile to a pickle

        As description.

        Args:
            path (str): The path to write the pickle to. (default:
            "configuration.params"`)
        """
        params = self.as_dict()

        with open(path, 'wb') as handle:
            pickle.dump(params, handle)

    def explain(self, param: str):
        """Disambiguates a parameter.

        Given a parameter, if possible, return a human-readable description
        of the parameter, including options. This is intended as an interactive
        diagnostic for any runs of kaitiaki which are interactive or are
        otherwise non-clustered.

        Args:
            param (str): The parameter to disambiguate.

        Notes:
            Writes to standard output, either through the debugger (if the
            parameter cannot be disambiguated) or through the print statement
            (if the parameter *can* be disambiguated).
        """
        val = self.get(param)
        param = param.lower()

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
        else:
            kaitiaki.debug('info',
                           f'Parameter {param} not disambiguable. Check the manual.')

    def show_in_file(self, param):
        """Shows a highlighted value in the datafile at stdout.

        Reprints faithfully the contents of the datafile, with a given
        parameter highlighted in bright green (as per colorama).

        Args:
            param (str): The parameter to show

        Raises:
            KeyError: If the parameter does not exist in the datafile

        Notes:
            Prints a rather large amount of text to stdout, including text
            highlighted with colorama (in the form
            f"{Fore.GREEN}{Style.BRIGHT}{string}{Style.RESET_ALL}").
        """
        param = param.lower()
        if param not in kaitiaki.constants.dfile_struct.keys():
            raise KeyError('Invalid parameter')

        with self as data:
            # Now we're cooking
            colorama.init()
            # Find where in the datafile the parameter is:
            row_loc, start, finish = kaitiaki.constants.dfile_struct[param]
            substr = data._contents[row_loc][start:finish]

            # substr is now the contents of the datafile at the parameter
            # location (i.e., the value we want to highlight)

            infos = kaitiaki.constants.dfile_struct.values()
            j = max([info[0] for info in infos])
            # j is the maximum row number we need to print to to represent
            # the contents of the datafile (as we do not need to print
            # commentary that may come after the "important" parts of the
            # file)

            def hl(string):
                """Highlights a string using colorama.

                Args:
                    string (str): The string to highlight

                Returns:
                    str: The highlighted string, (in bright green)
                """
                return Fore.GREEN + Style.BRIGHT + string + Style.RESET_ALL

            for i, row in enumerate(data._contents):
                # for each line...
                if i > j:
                    # ... are we past the meaningful bit? If so, stop.
                    break

                if i == row_loc:
                    # ...does this row contain the string we are looking for?
                    if finish is None:
                        # ... and is the string at the end of the line?
                        end = ''
                    else:
                        end = row[finish:]

                    # ... if it is in this row, print the highlighted string
                    print(row[:start] + hl(substr) + end)
                else:
                    # ...otherwise, just print the regular row.
                    print(row)

    def as_dict(self):
        """Fetches the entire datafile as a dictionary

        Returns a dictionary of {key: value} pairs representing the entries
        in the datafile.

        Returns:
            dict: THe dictionary containing the entries in the datafile.
        """
        keys = list(kaitiaki.constants.dfile_struct.keys())

        return self.get(keys)

    def __eq__(self, other):
        my_values = self.as_dict()
        other_val = other.as_dict()

        return my_values == other_val

    def changes_from_base(self):
        """Determine the changes between the current instance of the datafile
        and the default one.

        Runs a simple comparator of the present datafile against the default
        one that kaitiaki ships with.

        Returns:
            dict: A dictionary containing the mismatches for a given parameter.

        Notes:
            Consider a datafile where the only changes from base are that
            PARAM_A is set to 7 instead of 3, and PARAM_B is set to 3.4 instead
            of 3.0. Then the returned dictionary looks like:

            >>> {
            >>>     'PARAM_A': {'current': 7, 'base': 3}
            >>>     'PARAM_B': {'current': 3.4, 'base': 3.0}
            >>> }
        """
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
        """Returns the original string representation of the data file"""
        with self as data:
            return '\n'.join(data._contents)

    def compare(self, other, tag_names=('other', 'self')):
        """Compare two datafiles against each other.

        Element-wise comparison between two data files. The datafiles will be
        returned in a format (assuming default parameters, and that PARAM_A is
        4 in self and 5 in other):

        >>> {'PARAM_A': 'self': 4, 'other': 5}

        Args:
            other (DataFileParser): The comparative datafile
            tag_names (tuple): The names to assign the datafiles
                               (default: `('other', 'self')`)

        Returns:
            dict: The mismatches in the datafiles.
        """
        mismatches = dict()

        for key in kaitiaki.constants.dfile_struct.keys():
            other_val = other.get(key)
            my_val = self.get(key)

            if other_val != my_val:
                mismatches[key] = {tag_names[0]: other_val,
                                   tag_names[1]: my_val}

        return mismatches

    def __enter__(self):
        """Entry method for the DataFileParser as a context manager."""
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
                """Configures a datafile from a given pickle.

                Given a pickled file, set all of the parameters in the pickle
                as parameters in the datafile. Designed to be an interoperable
                way of translating datafiles across projects.

                Args:
                    path (str): The path to the pickle to load.
                """
                with open(path, 'rb') as handle:
                    params = pickle.load(handle)

                self.set(params)

            def _check_scientific_notation(self, param):
                """Utility method: Compute the number of s.f. required for
                scientific notation of a given datafile parameter.

                Returns the number of significant figures required for a given
                parameter  in the datafile, that is represented as scientific
                notation. The number of significant figures in the number
                represented as mEn, i.e., m*10^n, is m, and we return m-1 --
                that is, we do not include the digit to the left of the
                decimal point.

                Args:
                    param (str): The parameter to lookup

                Returns:
                    number: The number of significant figures.
                """
                idx = self._get_index_of_parameter(param)
                if idx[0] == 3:
                    # This line requires scientific notation to 1sf
                    return 1
                elif idx[0] in [17, 18]:
                    # These lines require scientific notation to 2sf
                    return 2
                else:
                    param = param.lower()
                    if param in ['ct8', 'ct9', 'ct10']:
                        return 1
                    if param == 'zs':
                        return 2
                    if param in ['vrot1', 'vrot2']:
                        return 2
                    if param in ['facsgmin', 'sgthfac']:
                        return 2
                    if param in ['hkh', 'gff']:
                        return 2

                return 0

            def _determine_decimal_places(self, param):
                """Utility method: Determine the number of decimal places
                for a given parameter NOT in scientific notation.

                Determines the number of decimal places (to the RIGHT of the
                decimal point) for a given parameter.

                Args:
                    param (str): The parameter to lookup

                Returns:
                    number: The number of decimal places.
                """
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
                """Utility function: Returns the location in the datafile of a
                parameter.

                Returns a 3-tuple of the location in the datafile corresponding
                to a given parameter. The 3-tuple is of the form
                (line_number, start_position, end_position).

                Args:
                    param (str): The parameter to lookup

                Returns:
                    tuple(int, int, int): The location of the parameter

                Raises:
                    KeyError: If an invalid parameter is passed.
                """

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
                """Utility function: See self.elucidate()."""
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

                """Generate a report of all of the parameters in the data file
                as key-value pairs.

                Creates a file with key: value pairs of parameters. This file
                is intended to be human-readable rather than machine-readable
                (you should be using the raw datafile for machine-readable
                methods anyways).

                Args:
                    file (str): The file to output the elucidation to.

                Keyword Args:
                    eludicate_only (list): A list of parameters to elucidate.
                                           If None, then all parameters in the
                                           datafile will be elucidated.
                                           (default: None)

                    fmt (str): The file format to elucidate (either `tex` or
                               `plaintext`). `tex` creates a booktabs table
                               to be pasted into a TeX document, whereas
                               `plaintext` is just a list of key-value pairs.

                    splice_into (int): An integer representing how many
                                       column-pairs to splice the list into if
                                       we are printing a tex'd table. For
                                       instance, splice_int0=2 gives a table
                                       like: Parameter | Value | Comment | Parameter | Value | Comment |
                """
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

            def _write_to_pointer(self, pointer, data, mode: str = 'w'):
                """Writes to a file handler object.

                Writes to a file handler object (called a pointer here for
                backwards compatibility reasons).

                Args:
                    pointer (file): the file object to write to (e.g., call
                                    self._write_to_pointer(open('some_file',
                                    'w'), ...) or similar).
                    data (str): The data to write to the pointer.
                    mode (str): Whether to "replace" or append ("w") the file.
                """
                options = ['w', 'replace']
                error = "Invalid write mode selected. Valid options are: "
                error += ", ".join(options)

                assert mode in options, error

                if mode == 'replace':
                    pointer.seek(0)

                pointer.write("\n".join(data))

            def make_backup(self):
                """Backs up the datafile.

                Will always overwrite the contents of the backup datafile. The
                backup datafile is assumed to be the file "[file].bak", where
                [file] is the name of the datafile, in the same directory as
                the datafile. For the average user, this is the file data.bak
                in the run directory.

                As a general rule, since this forcibly overwrites the data in
                data.bak, you shouldn't use this method; you should use
                self.backup_if_not_exists() instead (which *doesn't* overwrite)
                """
                with open(self._datafile + ".bak", 'w') as file:
                    self._write_to_pointer(file, self._original_contents, 'w')

            def backup_if_not_exists(self):
                """Creates a backup of the datafile if one doesn't already
                exist.

                Companion method to self.make_backup(). This one checks first
                to see if the targeted backup file exists or not, and only
                creates a backup if one doesn't already exist.
                """
                if not path.exists(self._datafile + ".bak"):
                    self.make_backup()

            def _setitem(self, param, value):
                """Sets the parameter param to value value in the datafile.

                Writes a given parameter to the datafile. This method makes
                changes to disk; it commits the parameter and doesn't just
                remember it. Thus, a parameter set by this method is available
                to the STARS code itself.

                Args:
                    param (str): The parameter to set the value of
                    value (mixed): The value to set the parameter to.
                """
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

                # Format the value into what the STARS code expects.
                if scientific_notation_dp > 0:
                    value = f'%.{scientific_notation_dp}E' % Decimal(value)
                if num_dp > 0:
                    value = f'%.{num_dp}f' % Decimal(value)

                # Determine where in the datafile to stop writing to
                if idx[2] is None:
                    # The end of the line?
                    endpoint = len(self._contents[idx[0]])
                else:
                    # Or where we are told the endpoint is?
                    endpoint = idx[2]

                # Partition this line into
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
