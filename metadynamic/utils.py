from multiprocessing import current_process
from time import process_time
from typing import Callable, TypeVar, Any, Optional
from logging import getLogger, FileHandler, StreamHandler, Handler, Logger
from datetime import datetime

T = TypeVar("T")
A = TypeVar("A")


def memoize_property(f: Callable[[T], Any]) -> Callable[[T], Any]:
    def memoized(self: T) -> Any:
        try:
            return self.__getattribute__(f"_{f.__name__}_memoized")
        except AttributeError:
            self.__setattr__(f"_{f.__name__}_memoized", f(self))
            return self.__getattribute__(f"_{f.__name__}_memoized")

    return memoized


def memoize_oneparam(f: Callable[[T, A], Any]) -> Callable[[T, A], Any]:
    def memoized(self: T, pos: A) -> Any:
        try:
            return self.__getattribute__(f"_{f.__name__}_{pos}_memoized")
        except AttributeError:
            self.__setattr__(f"_{f.__name__}_{pos}_memoized", f(self, pos))
            return self.__getattribute__(f"_{f.__name__}_{pos}_memoized")

    return memoized


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def get_all_subclasses(cls):  # Annotation ?????
    all_subclasses = []
    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))
    return all_subclasses


class Timer:
    def __init__(self):
        self.reset()

    def reset(self):
        self._ptime0 = process_time()

    @property
    def time(self) -> float:
        return process_time() - self._ptime0

    
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


class HappyEnding(Finished):
    pass


class BadEnding(Finished):
    pass


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


class RuntimeLim(Finished):
    num = 5
    error_message = "Runtime limit exceeded"


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
    def time():
        return datetime.now().strftime("%H:%M:%S, %d/%m/%y")

    def __init__(
        self, timer: Timer, filename: Optional[str] = None, level: str = "INFO"
    ):
        if filename:
            if filename.count(".") != 1:
                raise ValueError("Please enter filename as 'filename.log'")
            basename, suf = filename.split(".")
            self.filenamethread = basename + "-{}." + suf
        self._timer = timer
        self.filename = filename
        self.level = level
        self.connected = False
        self.connect("Logger creation")

    def connect(self, reason: str = "unknown", thread: Optional[int] = None):
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

    def disconnect(self, reason: str = "unknown"):
        if self.connected:
            self.debug(f"Disconnecting; reason: {reason}")
            assert self._handler is not None
            self._logger.removeHandler(self._handler)
            self._handler.close()
            self._logger = BlackholeLogger("BlackHole")
            self._handler = None
            self.connected = False

    def debug(self, msg: str):
        self._logger.debug(
            f"DEBUG-{current_process().name} : {msg}   (rt={self._timer.time}, t={self.time()})"
        )

    def info(self, msg: str):
        self._logger.info(
            f"INFO-{current_process().name} : {msg}   (rt={self._timer.time}, t={self.time()})"
        )

    def warning(self, msg: str):
        self._logger.warning(
            f"WARNING-{current_process().name} : {msg}   (rt={self._timer.time}, t={self.time()})"
        )

    def error(self, msg: str):
        self._logger.error(
            f"ERROR-{current_process().name} : {msg}   (rt={self._timer.time}, t={self.time()})"
        )

