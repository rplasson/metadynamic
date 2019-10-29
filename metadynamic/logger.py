from multiprocessing import current_process
from time import process_time
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime
from typing import Optional


class Timer:
    def __init__(self):
        self.reset()

    def reset(self):
        self._ptime0 = process_time()

    @property
    def time(self) -> float:
        return process_time() - self._ptime0


class BlackholeLogger(Logger):
    def debug(self, msg: str, *args, **kwargs):
        pass

    def info(self, msg: str, *args, **kwargs):
        pass

    def warning(self, msg: str, *args, **kwargs):
        pass

    def error(self, msg: str, *args, **kwargs):
        pass


class Log:
    @staticmethod
    def time() -> str:
        return datetime.now().strftime("%H:%M:%S, %d/%m/%y")

    def __init__(self, filename: Optional[str] = None, level: str = "INFO"):
        if filename:
            if filename.count(".") != 1:
                raise ValueError("Please enter filename as 'filename.log'")
            basename, suf = filename.split(".")
            self.filenamethread: str = basename + "-{}." + suf
        self._timer: Timer = Timer()
        self.filename: str = filename
        self.level: str = level
        self.connected: bool = False
        self.connect("Logger creation")

    def connect(self, reason: str = "unknown", thread: Optional[int] = None) -> None:
        if not self.connected:
            filename = (
                self.filenamethread.format(thread)
                if thread and self.filename
                else self.filename
            )
            self._logger = getLogger("Polymer Log")
            self._handler: Optional[Handler] = FileHandler(
                filename
            ) if filename else StreamHandler()
            self._logger.addHandler(self._handler)
            self._logger.setLevel(self.level)
            self.debug(f"Connected to {filename}; reason: {reason}")
            self.connected = True

    def disconnect(self, reason: str = "unknown") -> None:
        if self.connected:
            self.debug(f"Disconnecting; reason: {reason}")
            assert self._handler is not None
            self._logger.removeHandler(self._handler)
            self._handler.close()
            self._logger = BlackholeLogger("BlackHole")
            self._handler = None
            self.connected = False

    def debug(self, msg: str) -> None:
        self._logger.debug(
            f"DEBUG-{current_process().name} : {msg}   (rt={self.runtime()}, t={self.time()})"
        )

    def info(self, msg: str) -> None:
        self._logger.info(
            f"INFO-{current_process().name} : {msg}   (rt={self.runtime()}, t={self.time()})"
        )

    def warning(self, msg: str) -> None:
        self._logger.warning(
            f"WARNING-{current_process().name} : {msg}   (rt={self.runtime()}, t={self.time()})"
        )

    def error(self, msg: str) -> None:
        self._logger.error(
            f"ERROR-{current_process().name} : {msg}   (rt={self.runtime()}, t={self.time()})"
        )

    def runtime(self) -> float:
        return self._timer.time

    def reset_timer(self) -> None:
        self._timer.reset()
