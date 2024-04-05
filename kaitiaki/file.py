import kaitiaki

import os
from .file_handlers import plotfile, outfile, datafile, modelfile, sneplotfile

plot = plotfile.plot
plot2 = plotfile.plot2
sneplot = sneplotfile.sneplot
sneplot2 = sneplotfile.sneplot2

out = outfile.outfile
out2 = outfile.outfile2
outfile = outfile

data = datafile.DataFileParser

# Syntactic sugar
modin = modelfile.modin
modout = modelfile.modout
modin2 = modelfile.modin2
modout2 = modelfile.modout2


def load(path_to_file: str):
    """Loads a STARS file.

    Loads a STARS file (that has not been renamed).

    Args:
        path_to_file (str): The path to the file to load

    Returns:
        (mixed): An instance of the file

    Raises:
        FileNotFoundError: If the file you are attempting to load doesn't exist
        ValueError: If the file type cannot be inferred from the filename.
    """
    if not os.path.exists(path_to_file):
        raise FileNotFoundError(f"{path_to_file} doesn't exist.")

    filename_to_load = path_to_file.split(os.sep)[-1]

    match filename_to_load:
        case 'plot' | 'plot2':
            return plot(path_to_file)
        case 'out' | 'out2':
            return out(path_to_file)
        case 'modin' | 'modin2':
            return modin(path_to_file)
        case 'modout' | 'modout2':
            return modout(path_to_file)
        case 'data':
            kaitiaki.debug('info', ('Loading the data file returns a '
                                    'read-only object. You cannot '
                                    'change datafile parameters with '
                                    'this.'))
            return data(path_to_file)
        case _:
            raise ValueError(('Could not figure out what file this is '
                              'supposed to be. Please use one of the '
                              'direct initialisation methods.'))
