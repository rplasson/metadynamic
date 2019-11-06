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


class InputError(Finished):
    pass


class FileNotFound(InputError):
    num = 6
    error_message = "The provided file was not found."


class BadFile(InputError):
    num = 7
    error_message = "The provided file is badly formed"
