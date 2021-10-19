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

from kaitiaki._metadata import __version__

def debug(msgtype, message, fatal=True):
    """General purpose debug message handler

    Allows us to print to stdout when debugging (developing) and fail
    on production.

    Arguments:
        msgtype {string} -- The message type to throw. Must be 'info', 'warning', or 'error'.

        message {string} -- The message to throw.

    Keyword Arguments:
        fatal {bool} -- Whether or not the message should be a fatal
                        error. Ignored if takahe.constants.DEBUG_MODE
                        is True. (default: {True})

    Raises:
        takahe.TakaheWarning    -- A warning type if we are not in debug
                                   mode and the error should be fatal.
        takahe.TakaheFatalError -- An error type if we are not in debug
                                   mode and the error should be fatal.
    """
    assert msgtype in ['warning', 'error', 'info'], ("Message type "
                                                     f"\"{msgtype}\" is "
                                                     "not recognised.")

    assert fatal == True or fatal == False, "fatal must be a boolean type."

    if constants.DEBUG_MODE:
        if msgtype == 'warning':
            header = "\033[93m\033[1m[WARNING]\033[0m "
        elif msgtype == 'error':
            header = "\033[91m\033[1m[ERROR]\033[0m "
        elif msgtype == 'info':
            header = "\033[96m\033[1m[INFO]\033[0m "

        print(header + str(message))
