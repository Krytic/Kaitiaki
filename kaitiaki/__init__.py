#
#
#      .                 .--.         .--.              .     .--.
#    .'|                 |__|         |__|            .'|     |__|
#  .'  |                 .--.     .|  .--.          .'  |     .--.
# <    |            __   |  |   .' |_ |  |    __   <    |     |  |
#  |   | ____    .:--.'. |  | .'     ||  | .:--.'.  |   | ____|  |
#  |   | \ .'   / |   \ ||  |'--.  .-'|  |/ |   \ | |   | \ .'|  |
#  |   |/  .    `" __ | ||  |   |  |  |  |`" __ | | |   |/  . |  |
#  |    /\  \    .'.''| ||__|   |  |  |__| .'.''| | |    /\  \|__|
#  |   |  \  \  / /   | |_      |  '.'    / /   | |_|   |  \  \
#  '    \  \  \ \ \._,\ '/      |   /     \ \._,\ '/'    \  \  \
# '------'  '---'`--'  `"       `'-'       `--'  `"'------'  '---'
#
#  _                                         _      _                   _
# | |                                       (_)    | |                 | |
# | |__  _   _   ___  ___  __ _ _ __    _ __ _  ___| |__   __ _ _ __ __| |___
# | '_ \| | | | / __|/ _ \/ _` | '_ \  | '__| |/ __| '_ \ / _` | '__/ _` / __|
# | |_) | |_| | \__ \  __/ (_| | | | | | |  | | (__| | | | (_| | | | (_| \__ \
# |_.__/ \__, | |___/\___|\__,_|_| |_| |_|  |_|\___|_| |_|\__,_|_|  \__,_|___/
#         __/ |
#        |___/


#                   .
#
#                    .
#          /^\     .
#     /\   "V"
#    /__\   I      O  o
#   //..\\  I     .
#   \].`[/  I
#   /l\/j\  (]    .  O
#  /. ~~ ,\/I          .
#  \\L__j^\/I       o
#   \/--v}  I     o   .
#   |    |  I   _________
#   |    |  I c(`       ')o
#   |    l  I   \.     ,/
# _/j  L l\_!  _//^---^\\_
#
#      Here be wizard.

# Core Python
import os
import importlib.resources as importlib_resources

# Core kaitiaki
import kaitiaki.STARSController as STARS
import kaitiaki.helpers as helpers
import kaitiaki.constants as constants
import kaitiaki.classifier as classify
import kaitiaki.OptionLexer as lexer
import kaitiaki.kicks as kicks
import kaitiaki.sntools as sntools

# Augments -- will probably be made its own project
import kaitiaki.augments as augments

# Utils
import kaitiaki.terminal as terminal
import kaitiaki.file as file
from .utils import transforms

# Other tools I wrote
import glisten

# Metadata
from kaitiaki._metadata import __version__

# Sugar
stars = STARS.STARSController


def load_file(filename):
    # data_path = pkg_resources.resource_filename('kaitiaki', '../backup_data')

    # pathname = os.path.join(data_path, '..', 'backup_data')
    # pathname = os.path.normpath(pathname)
    pathname = importlib_resources.files('kaitiaki') / '..' / 'backup_data'
    pathname = os.path.normpath(pathname)

    if '..' in filename:
        raise ValueError("Traversing up the filetree is not permitted.")

    allowed_folders = ['data.bak', 'COtables', 'dat', 'modins']

    for folder in allowed_folders:
        if filename.startswith(folder):
            break
    else:
        error = f"Folder not valid (allowed: {'/'.join(allowed_folders)})"
        raise ValueError(error)

    if not os.path.exists(os.path.join(pathname, filename)):
        raise IOError(f"{filename} is not a valid kaitiaki internal file.")

    with open(os.path.join(pathname, filename), 'rb') as f:
        file = f.read().decode('utf-8')

    return file


def format_metallicity(Z):
    if isinstance(Z, float):
        if Z in [1e-4, 1e-5]:
            return {
                1e-4: 'zem4',
                1e-5: 'zem5'
            }[Z]
        else:
            Z = str(Z).split('.')[1]
            Z = 'z' + Z.ljust(3, '0')
            return Z
    if isinstance(Z, str):
        if Z[0] != 'z':
            return format_metallicity(float(Z))
        return Z
    raise ValueError(f"I don't know what {Z} means.")


log = glisten.log.Logger('~/.kaitiaki-log', line_length_break=72*3)


def debug(msgtype, message, fatal=False):
    """General purpose debug message handler

    Allows us to print to stdout when debugging (developing) and fail
    on production.

    Arguments:
        msgtype {string} -- The message type to throw. Must be 'info', 'warning', or 'error'.

        message {string} -- The message to throw.

    """

    assert_err = f"Message type \"{msgtype}\" is not recognised."

    assert msgtype in ['warning',
                       'error',
                       'info',
                       'status'], assert_err

    match msgtype:
        case 'info':
            log.info(message)
        case 'status':
            log.status(message)
        case 'warning':
            log.warn(message)
        case 'error':
            log.error(message)


def deprecate(func):
    return log.deprecate(func)


def run(lexer):
    STARS = kaitiaki.STARS.STARSController()

    for options in lexer:
        STARS.configure_parameters(options)
