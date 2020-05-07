#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2019 by Raphaël Plasson
#
# This file is part of metadynamic
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

"""General exceptions to be used in metadynamic code.

Exception categories inherit from Finished:
 - L{HappyEnding} is intended for normal end of the computation
 - L{Aborted} is intended to signal that the computation ended earlier than expected,
   but lead to a nonetheless correct result
 - L{BadEnding} is intended to signal that something when wrong during the computation
 - L{InputError} is intended to signal problems with file read or write, unrelated to
   the computation itself.

All exceptions intended to be used are derived from one of these categories, indicating the general
context of the exception.

"""

from typing import Union, Callable
from types import FrameType
from signal import (
    signal,
    getsignal,
    Signals,
    Handlers,
    SIGTERM,
    SIGINT,
    SIG_IGN,
    SIG_DFL,
)


class Finished(Exception):
    """General exception raised at run completion."""

    num = -1
    """exception number"""
    error_message = "Stopped for unknown reason"
    """exception base message"""

    def __init__(self, detail: str = ""):
        """When the exception is raised, details can be added to the base error mesage.

        @param detail: details to be added to error message
        @type detail: str

        """
        self.detail = detail
        """additional information"""
        super().__init__()

    @property
    def message(self) -> str:
        """Format the exception message.

        @return: full exception message
        @rtype: str

        """
        msg = self.error_message
        if self.detail != "":
            msg = msg + " -> " + self.detail
        return msg

    def __str__(self) -> str:
        """Return a formatted ending message."""
        return f"End ({self.num}): {self.message}"


# Exception categories


class HappyEnding(Finished):
    """Raised for normal run completion."""


class BadEnding(Finished):
    """Raised for faulty run completion."""


class Aborted(Finished):
    """Raised for shortened, but correct, run completion."""


class InputError(BadEnding):
    """Raised when problem are encountered when reading input files."""


# Exceptions expected to be raised

# Generic


class InternalError(BadEnding):
    """Something went bad in the code."""

    num = 0
    error_message = "Something went bad in the code"


# Regular ending


class TimesUp(HappyEnding):
    """Time is up."""

    num = 10
    error_message = "Time is up"


class NoMore(HappyEnding):
    """No more reactions can be processed."""

    num = 11
    error_message = "No more reactions can be processed"


# Problems during run


class NotFound(BadEnding):
    """No reaction could be find."""

    num = 20
    error_message = "No reaction could be find"


class RoundError(BadEnding):
    """Rounding problem/negative probability."""

    num = 21
    error_message = "Rounding problem/negative probability detected"


class DecrZero(BadEnding):
    """Tried to decrement unpopulated species."""

    num = 22
    error_message = "Tried to decrement unpopulated species"


# Shortened runs


class RuntimeLim(Aborted):
    """Runtime limit exceeded."""

    num = 30
    error_message = "Runtime limit exceeded"


class InitError(Aborted):
    """Error during initialization."""

    num = 31
    error_message = "Error during initialization"


class Interrupted(Aborted):
    """Asked to stop."""

    num = 32
    error_message = "Asked to stop"


class OOMError(Aborted):
    """Out of Memory."""

    num = 33
    error_message = "Out of Memory"


# IO errors


class FileNotFound(InputError):
    """The provided file was not found."""

    num = 40
    error_message = "The provided file was not found."


class BadFile(InputError):
    """The provided file is badly formed."""

    num = 41
    error_message = "The provided file is badly formed"


class BadJSON(InputError):
    """Bad JSON format."""

    num = 42
    error_message = "Bad JSON format"


class FileCreationError(InputError):
    """The file couldn't be created."""

    num = 43
    error_message = "The file couldn't be created"


class NotAFolder(InputError):
    """The provided foldername is not a folder."""

    num = 44
    error_message = "The provided foldername is not a folder"


# Signal handling

SignHandler = Union[Callable[[Signals, FrameType], None], int, Handlers, None]
"""Generic type that can be returned as a signal handler"""


class SignalCatcher:
    """Catch temporarily SIGTERM and SIGINT signal interruptions."""

    def __init__(self) -> None:
        """Simply create the object, and do not change anything to signal handling."""
        self.alive: bool = False
        """aliveness flag
        (when set in listen state, it is True until SIGINT or SIGTERM is received)"""
        self.signal: str = ""
        """Name of the signal received while in 'listen' state"""
        self.frame: str = ""
        """Frame of the signal received while in 'listen' state"""
        self._initial_term: SignHandler = getsignal(SIGTERM)
        """Initial handler to which SIGTERM was connected"""
        self._initial_int: SignHandler = getsignal(SIGINT)
        """Initial handler to which SIGTINT was connected"""

    def reset(self) -> None:
        """Return to the signal handling state as it was at object creation."""
        signal(SIGTERM, self._initial_term)
        signal(SIGINT, self._initial_int)

    def ignore(self) -> None:
        """SIGTERM and SIGINT signals to be ignored."""
        self.init_signal(SIG_IGN)

    def release(self) -> None:
        """SIGTERM and SIGINT signal to be set to default behaviour."""
        self.init_signal(SIG_DFL)

    def listen(self) -> None:
        """SIGTERM and SIGINT signal to be set to listen.

        When a signal is received, L{signal_listen} is called, namely setting the flag alive to
        False

        """
        self.alive = True
        self.init_signal(self.signal_listen)

    @staticmethod
    def init_signal(handler: SignHandler) -> None:
        """Connect the SIGTERM and SIGINT signals to a specific handler.

        @param handler: signal handler to connect to
        @type handler: SignHandler

        """
        signal(SIGTERM, handler)
        signal(SIGINT, handler)

    def signal_listen(self, received_signal: Signals, frame: FrameType) -> None:
        """Actions to be performed when SIGINT or SIGTERM are received when in 'listen' state.

        It switches the self.alive flag to False,
        saves the signal name in self.signal,
        saves the signal frame in self.frame,
        then switches the signal reception to 'ignore'

        @param received_signal: received signal
        @type received_signal: Signals
        @param frame: context frame at signal reception
        @type frame: FrameType

        """
        self.alive = False
        self.signal = Signals(received_signal).name
        self.frame = str(frame)
        self.ignore()
