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

"""
metadynamic.logger
==================

Logging facility, taking into account logging from a MPI run, and the
possibility to save logs in a hdf5 file simultaneously to normal logging.


Provides
--------

 - L{Timer}: Simple class for tracking passed processing time from a checkpoint

 - L{Log}: High level class for MPI-aware logging, with log save in hdf5.

 - L{LOGGER}: Global Log object

"""

from time import process_time
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime
from typing import Optional

from metadynamic.inval import invalidstr, isvalid
from metadynamic.hdf5 import ResultWriter
from metadynamic.mpi import MPI_STATUS


class Timer:
    """Simple class for tracking passed processing time from a checkpoint"""

    def __init__(self) -> None:
        self._ptime0: float
        """Checkpoint time"""
        self.reset()

    def reset(self) -> None:
        """Set the checkpoint at the present time"""
        self._ptime0 = process_time()

    @property
    def time(self) -> float:
        """Passed processing time from the checkpoint"""
        return process_time() - self._ptime0


class Log:
    """High level class for MPI-aware logging, with log save in hdf5."""

    def __init__(
        self,
        filename: str = invalidstr,
        level: str = "INFO",
        timeformat: str = "%H:%M:%S, %d/%m/%y",
    ):
        """
        @param filename: name of log file
        @type filename: str
        @param level: logging level (INFO, DEBUG, etc.) (Default value = "INFO")
        @type level: str
        @param timeformat: format string for date/time output (Default value = "%H:%M:%S, %d/%m/%y")
        @type timeformat: str
        """
        self.connected: bool = False
        """Connected flag (is the logger connected to a Handler?)"""
        self.timeformat: str = timeformat
        """format string for date/time output"""
        self._timer: Timer = Timer()
        self._logger: Logger = getLogger("Metadynamic Log")
        self._handler: Handler
        self.level: str
        """logging level"""
        self.filename: str
        """output log file"""
        self.writer: Optional[ResultWriter] = None
        """writer to hdf5 file"""
        self.setlevel(level)
        self.connect(filename)

    def setsaver(self, writer: ResultWriter) -> None:
        """Connect to a hdf5 writer

        @param writer: hdf5 writer
        @type writer: ResultWriter
        """
        self.writer = writer

    def setlevel(self, level: str = "INFO") -> None:
        """Set logging level

        @param level: logging level (Default value = "INFO")
        @type level: str
        """
        self.debug(f"Switched to level {level}")
        self.level = level
        self._logger.setLevel(self.level)

    def connect(self, filename: str = "") -> None:
        """direct logging to textfile 'filename'
        if empty, log to standard output

        @param filename: log file name
        @type filename: str
        """
        self.filename = filename if filename else invalidstr
        dest = filename if isvalid(filename) else "stream"
        if self.connected:
            self.disconnect("Disconnecting old handler before connection to new one.")
        self._handler = FileHandler(self.filename) if isvalid(self.filename) else StreamHandler()
        self._logger.addHandler(self._handler)
        self.debug(f"Connected to {self.filename}; reason: Logger directed to {dest}")
        self.connected = True

    def disconnect(self, reason: str = "unknown") -> None:
        """Disconnect to current handler

        @param reason: message to debug log (Default value = "unknown")
        @type reason: str
        """
        self.debug(f"Disconnecting; reason: {reason}")
        self._logger.removeHandler(self._handler)
        self._handler.close()
        self.connected = False

    def _format_msg(self, level: str, msg: str) -> str:
        """Format the log message.

        @param level: severity level
        @type level: str
        @param msg: log message
        @type msg: str
        @return: formatted message
        @rtype:str
        """
        return (
            f"{level}-{MPI_STATUS.rank} : {msg}   (rt={self.runtime}, t={self.time})"
        )

    def _savelog(self, level: int, msg: str) -> None:
        """Save the log message to hdf5.

        @param level: severity level
        @type level: int
        @param msg: log message
        @type msg: str
        """
        if self.writer is not None:
            self.writer.write_log(level, self.time, self.runtime, msg)

    def debug(self, msg: str) -> None:
        """log message at severity 'DEBUG'

        @param msg: log message
        @type msg: str
        """
        self._savelog(10, msg)
        self._logger.debug(self._format_msg("DEBUG", msg))

    def info(self, msg: str) -> None:
        """log message at severity 'INFO'

        @param msg: log message
        @type msg: str
        """
        self._savelog(20, msg)
        self._logger.info(self._format_msg("INFO", msg))

    def warning(self, msg: str) -> None:
        """log message at severity 'WARNING'

        @param msg: log message
        @type msg: str
        """
        self._savelog(30, msg)
        self._logger.warning(self._format_msg("WARNING", msg))

    def error(self, msg: str) -> None:
        """log message at severity 'ERROR'

        @param msg: log message
        @type msg: str
        """
        self._savelog(40, msg)
        self._logger.error(self._format_msg("ERROR", msg))

    @property
    def time(self) -> str:
        """Current time

        @rtype time: str
        """
        return datetime.now().strftime(self.timeformat)

    @property
    def runtime(self) -> float:
        """run time (in seconds)

        @rtype time: float
        """
        return self._timer.time

    def reset_timer(self) -> None:
        """reset run timer to 0"""
        self.debug("Will reset the timer")
        self._timer.reset()


LOGGER = Log()
"""global object for logging messages"""
