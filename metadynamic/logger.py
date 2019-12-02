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

from multiprocessing import current_process
from time import process_time
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime
from typing import Union

from metadynamic.inval import InvalidStr, InvalidInt, Invalid, isvalid


class InvalidHandler(Invalid, Handler):
    _invalrepr = "Invalid Handler"


class Timer:
    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self._ptime0 = process_time()

    @property
    def time(self) -> float:
        return process_time() - self._ptime0


class DummyTimer(Timer):
    def reset(self) -> None:
        self._ptime0 = 0.0

    @property
    def time(self) -> float:
        return 0


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


class Log:
    @staticmethod
    def time() -> str:
        return datetime.now().strftime("%H:%M:%S, %d/%m/%y")

    def __init__(self, filename: str = InvalidStr(), level: str = "INFO"):
        if isvalid(filename):
            if filename.count(".") != 1:
                raise ValueError("Please enter filename as 'filename.log'")
            basename, suf = filename.split(".")
            self.filenamethread: str = basename + "-{}." + suf
        self.filename: str = filename
        self._timer: Timer
        self._logger: Union[Logger, BlackholeLogger]
        self._handler: Handler
        self.level: str = level
        self.connected: bool = False
        self.connect("Logger creation")

    def connect(self, reason: str = "unknown", thread: int = InvalidInt()) -> None:
        if not self.connected:
            self._timer = Timer()
            filename = (
                self.filenamethread.format(thread)
                if isvalid(thread) and isvalid(self.filename)
                else self.filename
            )
            self._logger = getLogger("Polymer Log")
            self._handler = (
                FileHandler(filename) if isvalid(filename) else StreamHandler()
            )
            self._logger.addHandler(self._handler)
            self._logger.setLevel(self.level)
            self.debug(f"Connected to {filename}; reason: {reason}")
            self.connected = True
        else:
            self.reset_timer()
            self.debug(f"Attempted to reconnect; reason: {reason}")

    def disconnect(self, reason: str = "unknown") -> None:
        if self.connected:
            self.debug(f"Disconnecting; reason: {reason}")
            self._timer = DummyTimer()
            self._logger.removeHandler(self._handler)
            self._handler.close()
            self._logger = BlackholeLogger()
            self._handler = InvalidHandler()
            self.connected = False
        else:
            self.debug(f"Attempted to redisconnect; reason: {reason}")

    def _format_msg(self, origin: str, msg: str) -> str:
        return f"{origin}-{current_process().name} : {msg}   (rt={self.runtime()}, t={self.time()})"

    def debug(self, msg: str) -> None:
        self._logger.debug(self._format_msg("DEBUG", msg))

    def info(self, msg: str) -> None:
        self._logger.info(self._format_msg("INFO", msg))

    def warning(self, msg: str) -> None:
        self._logger.warning(self._format_msg("WARNING", msg))

    def error(self, msg: str) -> None:
        self._logger.error(self._format_msg("ERROR", msg))

    def runtime(self) -> float:
        return self._timer.time

    def reset_timer(self) -> None:
        self._timer.reset()


class Logged:
    log: Log

    @classmethod
    def setlogger(cls, filename: str = InvalidStr(), level: str = "INFO") -> None:
        cls.log = Log(filename, level)
