from multiprocessing import current_process
from time import process_time
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime
from typing import Optional, Union


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

    def removeHandler(self, handler: Optional[Handler]) -> None:
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
        self.filename: Optional[str] = filename
        self._timer: Timer
        self._logger: Union[Logger, BlackholeLogger]
        self._handler: Optional[Handler]
        self.level: str = level
        self.connected: bool = False
        self.connect("Logger creation")

    def connect(self, reason: str = "unknown", thread: Optional[int] = None) -> None:
        if not self.connected:
            self._timer = Timer()
            filename = (
                self.filenamethread.format(thread)
                if thread and self.filename
                else self.filename
            )
            self._logger = getLogger("Polymer Log")
            self._handler = FileHandler(filename) if filename else StreamHandler()
            self._logger.addHandler(self._handler)
            self._logger.setLevel(self.level)
            self.debug(f"Connected to {filename}; reason: {reason}")
            self.connected = True

    def disconnect(self, reason: str = "unknown") -> None:
        if self.connected:
            self.debug(f"Disconnecting; reason: {reason}")
            self._timer = DummyTimer()
            assert self._handler is not None
            self._logger.removeHandler(self._handler)
            self._handler.close()
            self._logger = BlackholeLogger()
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
