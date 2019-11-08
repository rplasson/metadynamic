from typing import Callable, TypeVar, Any

T = TypeVar("T")
A = TypeVar("A")


def memoize_property(func: Callable[[T], Any]) -> Callable[[T], Any]:
    def memoized(self: T) -> Any:
        try:
            return self.__getattribute__(f"_{func.__name__}_memoized")
        except AttributeError:
            self.__setattr__(f"_{func.__name__}_memoized", func(self))
            return self.__getattribute__(f"_{func.__name__}_memoized")

    return memoized


def memoize_oneparam(func: Callable[[T, A], Any]) -> Callable[[T, A], Any]:
    def memoized(self: T, pos: A) -> Any:
        try:
            return self.__getattribute__(f"_{func.__name__}_{pos}_memoized")
        except AttributeError:
            self.__setattr__(f"_{func.__name__}_{pos}_memoized", func(self, pos))
            return self.__getattribute__(f"_{func.__name__}_{pos}_memoized")

    return memoized
