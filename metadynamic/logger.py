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

from os import path
from time import process_time
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime
from typing import Union, Optional

from metadynamic.inval import invalidstr, invalidint, Invalid, isvalid
from metadynamic.hdf5 import ResultWriter
from metadynamic.mpi import MPI_STATUS


class InvalidHandler(Invalid, Handler):
    _invalrepr = "Invalid Handler"


invalidhandler = InvalidHandler()


class Timer:
    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self._ptime0 = process_time()

    @property
    def time(self) -> float:
        return process_time() - self._ptime0


class BlackholeLogger:
    """a fake Logger object that does nothing"""

    def debug(self, msg: str) -> None:
        pass

    def info(self, msg: str) -> None:
        pass

    def warning(self, msg: str) -> None:
        pass

    def error(self, msg: str) -> None:
        pass

    def removeHandler(self, handler: Handler) -> None:
        pass

    def setLevel(self, level: str) -> None:
        pass


class Log:
    @staticmethod
    def time() -> str:
        return datetime.now().strftime("%H:%M:%S, %d/%m/%y")

    def __init__(self, filename: str = invalidstr, level: str = "INFO"):
        self._timer: Timer = Timer()
        self._logger: Union[Logger, BlackholeLogger] = BlackholeLogger()
        self._handler: Handler = invalidhandler
        self.level: str
        self.filename: str
        self.writer: Optional[ResultWriter] = None
        self.setlevel(level)
        self.settxt(filename)

    def setsaver(self, writer: ResultWriter) -> None:
        self.writer = writer

    def setlevel(self, level: str = "INFO") -> None:
        self.level = level
        self._logger.setLevel(self.level)

    def settxt(self, filename: str = invalidstr) -> None:
        if not filename:
            filename = invalidstr
        if isvalid(filename):
            if filename.count(".") != 1:
                raise ValueError("Please enter filename as 'filename.log'")
            basename, suf = path.splitext(filename)
            self.filenamethread: str = basename + "-{}" + suf
        self.filename = filename
        dest = filename if isvalid(filename) else "stream"
        self.connect(f"Logger directed to {dest}")

    def connect(self, reason: str = "unknown") -> None:
        self.disconnect("Reconnecting...")
        filename = (
            self.filenamethread.format(MPI_STATUS.rank)
            if isvalid(self.filename)
            else self.filename
        )
        self._logger = getLogger("Polymer Log")
        self._handler = (
            FileHandler(filename) if isvalid(filename) else StreamHandler()
        )
        self._logger.addHandler(self._handler)
        self._logger.setLevel(self.level)
        self.debug(f"Connected to {filename}; reason: {reason}")

    def disconnect(self, reason: str = "unknown") -> None:
        self.debug(f"Disconnecting; reason: {reason}")
        self._logger.removeHandler(self._handler)
        self._handler.close()
        self._logger = BlackholeLogger()
        self._handler = invalidhandler

    def _format_msg(self, origin: str, msg: str) -> str:
        return f"{origin}-{MPI_STATUS.rank} : {msg}   (rt={self.runtime()}, t={self.time()})"

    def savelog(self, level: int, msg: str) -> None:
        if self.writer is not None:
            self.writer.write_log(level, self.time(), self.runtime(), msg)

    def debug(self, msg: str) -> None:
        self.savelog(10, msg)
        self._logger.debug(self._format_msg("DEBUG", msg))

    def info(self, msg: str) -> None:
        self.savelog(20, msg)
        self._logger.info(self._format_msg("INFO", msg))

    def warning(self, msg: str) -> None:
        self.savelog(30, msg)
        self._logger.warning(self._format_msg("WARNING", msg))

    def error(self, msg: str) -> None:
        self.savelog(40, msg)
        self._logger.error(self._format_msg("ERROR", msg))

    def runtime(self) -> float:
        return self._timer.time

    def reset_timer(self) -> None:
        self._timer.reset()


LOGGER = Log()
