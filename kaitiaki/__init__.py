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

import kaitiaki.STARSController as STARS
import kaitiaki.GridMaker as gridmaker
import kaitiaki.helpers as helpers
import kaitiaki.constants as constants
import kaitiaki.classifier as classify
import kaitiaki.model as model
import kaitiaki.OptionLexer as lexer
import kaitiaki.kicks as kicks

import kaitiaki.terminal as terminal
import kaitiaki.file as file
import kaitiaki.file_handlers as _filehandler

import glisten

from kaitiaki._metadata import __version__

DEBUG_MODE = False

log = glisten.log.Logger('~/.kaitiaki-log', line_length_break=72*3)

def debug(msgtype, message, fatal=False):
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
    if msgtype == 'info':
        log.info(message)
    elif msgtype == 'status':
        log.status(message)
    elif msgtype == 'warning':
        log.warn(message)
    elif msgtype == 'error':
        log.error(message)

def deprecate(func):
    return log.deprecate(func)

def run(lexer):
    STARS = kaitiaki.STARS.STARSController()

    for options in lexer:
        STARS.configure_parameters(options)