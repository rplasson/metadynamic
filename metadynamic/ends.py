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

"""General exceptions to be used in metadynamic code

Exception categories inherit from Finished:
HappyEnding, BadEnding, Aborted, InputError

All exceptions intended to be used are derived from one
of these categories. indicating the general context of the
exception.
 - HappyEnding is intended for normal end of the computation
 - Aborted is intended to signal that the computation ended earlier
   than expected, but lead to a nonetheless correct result
 - HappyEnding is intended to signal that something when wrong during
   the computation
 - InputError is intended to signal problems with file read or write,
   unrekated to the computation itself.
"""

import signal


class Finished(Exception):
    num = -1
    error_message = "Stopped for unknown reason"

    def __init__(self, detail: str = ""):
        self.detail = detail
        super().__init__()

    @property
    def message(self) -> str:
        msg = self.error_message
        if self.detail != "":
            msg = msg + " -> " + self.detail
        return msg

    def __str__(self) -> str:
        return f"End ({self.num}): {self.message}"


# Exception categories


class HappyEnding(Finished):
    pass


class BadEnding(Finished):
    pass


class Aborted(Finished):
    pass


class InputError(Finished):
    pass


# Exceptions expected to be raised


class TimesUp(HappyEnding):
    num = 0
    error_message = "Time is up"


class NoMore(HappyEnding):
    num = 1
    error_message = "No more reactions can be processed"


class NotFound(BadEnding):
    num = 2
    error_message = "No reaction could be find"


class RoundError(BadEnding):
    num = 3
    error_message = "Rounding problem/negative probability detected"


class DecrZero(BadEnding):
    num = 4
    error_message = "Tried to decrement unpopulated species"


class RuntimeLim(Aborted):
    num = 5
    error_message = "Runtime limit exceeded"


class FileNotFound(InputError):
    num = 6
    error_message = "The provided file was not found."


class BadFile(InputError):
    num = 7
    error_message = "The provided file is badly formed"


class BadJSON(InputError):
    num = 8
    error_message = "Bad JSON format"


class InitError(Aborted):
    num = 9
    error_message = "Error during initialization"


class Interrupted(Aborted):
    num = 10
    error_message = "Asked to stop"


def signal_handler(received_signal, frame):
    raise Interrupted(f"Stopped by signal {received_signal} in frame {frame}")


def signal_ignore(received_signal, frame):
    pass


def init_signal(ignore: bool, listen=signal.SIGTERM):
    if ignore:
        signal.signal(listen, signal_ignore)
    else:
        signal.signal(listen, signal_handler)
