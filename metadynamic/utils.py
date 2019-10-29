from typing import Callable, TypeVar, Any

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
