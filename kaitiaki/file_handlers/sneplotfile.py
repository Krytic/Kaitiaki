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


class sneplot(kaitiaki.file_handlers.plotfile.plot):
    def __init__(self,
                 file: str = 'sneplot',
                 allow_pad_age: bool = True,
                 row: str = 'all',
                 dummy_object: bool = False):
        super().__init__(file, allow_pad_age, row, dummy_object)

    def parse_plotfile(self, fname: str = 'sneplot',
                       row: str = 'all',
                       is_dummy: bool = False):
        """Parses an sneplot file.

        Args:
            fname (str): The sneplot file to load
            row (str): The row to load (obsolute -- always pass "all")
            is_dummy (bool): Whether the dataframe should be empty or not

        Returns:
            [type]: [description]
        """

        c = kaitiaki.constants.SNEPLOT_FILE_COLUMNS

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
            # 100E16.9
            spec = ([2]
                    + [16 for _ in range(100)])
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
                                     widths=spec
                                     # infer_nrows=99999,
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


def sneplot2(file: str = 'sneplot2', *args, **kwargs):
    return sneplot(file, *args, **kwargs)
