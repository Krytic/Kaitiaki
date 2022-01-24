"""
                  .

                   .
         /^\     .
    /\   "V"
   /__\   I      O  o
  //..\\  I     .
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

import pkgutil

import kaitiaki.DataFileParser as datafile
import kaitiaki.plotfile as plotfile
import kaitiaki.STARSController as STARS
import kaitiaki.GridMaker as gridmaker

from kaitiaki._metadata import __version__

DEBUG_MODE = False

def debug(msgtype, message, logfile_location=".log", fatal=False):
    """General purpose debug message handler

    Allows us to print to stdout when debugging (developing) and fail
    on production.

    Arguments:
        msgtype {string} -- The message type to throw. Must be 'info', 'warning', or 'error'.

        message {string} -- The message to throw.

    """

    assert msgtype in ['warning',
                       'error',
                       'info',
                       'status'], f"Message type \"{msgtype}\" is not recognised."

    endc = "\033[0m"

    if msgtype == 'warning':
        hdr = "\033[93m\033[1m"
        header = "[WARNING]"
    elif msgtype == 'error':
        hdr = "\033[91m\033[1m"
        header = "[ERROR]"
    elif msgtype == 'info':
        hdr = "\033[96m\033[1m"
        header = "[INFO]"
    elif msgtype == 'status':
        hdr = "\033[92m\033[1m"
        header = "[PROGRAM STATUS]"

    if DEBUG_MODE or fatal:
        print(f"{hdr}{header}{endc} {message}")
    else:
        with open(logfile_location, 'a') as f:
            f.write(f"{header} {message}\n")