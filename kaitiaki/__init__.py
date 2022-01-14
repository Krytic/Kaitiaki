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

def debug(msgtype, message):
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

    if msgtype == 'warning':
        header = "\033[93m\033[1m[WARNING]\033[0m "
    elif msgtype == 'error':
        header = "\033[91m\033[1m[ERROR]\033[0m "
    elif msgtype == 'info':
        header = "\033[96m\033[1m[INFO]\033[0m "
    elif msgtype == 'status':
        header = "\033[92m\033[1m[PROGRAM STATUS]\033[0m "

    print(header + str(message))
