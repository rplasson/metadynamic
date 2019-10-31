from typing import Generic, TypeVar, Dict, Set, List

from metadynamic.logger import Logged

T = TypeVar("T")


class Collect(Generic[T], Logged):
    _colltype = "Generic"

    def __init__(self, system: "System", drop: bool = False, categorize: bool = True):
        self.system = system  # remove this dependency???
        self.pool: Dict[str, T] = {}
        self.categories: Dict[str, Set[T]] = {}
        self.active: Dict[str, T] = self.pool if drop else {}
        self.categorize = categorize

    def __repr__(self) -> str:
        return f"<Collect of {len(self.pool)} {self._colltype}>"

    def __getitem__(self, name: str) -> T:
        """ Return the object as described by its name
            If it is the first call of the object, create it
            Else, return the already created one"""
        try:
            return self.pool[name]
        except KeyError:
            newobj = self._create(name)
            self.pool[name] = newobj
            self._initialize(newobj)
            return newobj

    def activate(self, name: str) -> None:
        """Put the object 'name' in the active section, then categorize it"""
        obj = self[name]
        # will fail if activate an duplicated object (i.e. 2 different objects with same name exists)
        assert obj is self[name]
        if name not in self.active:
            self.active[name] = obj
        if self.categorize:
            for catname in self._categorize(obj):
                try:
                    self.categories[catname].add(obj)
                except KeyError:
                    self.categories[catname] = {obj}

    def unactivate(self, name: str) -> None:
        """Remove the object 'name' from the active section, then
           uncategorize it"""
        try:
            del self.active[name]
        except KeyError:
            pass
        if self.categorize:
            for cat in self.categories.values():
                try:
                    cat.remove(self[name])
                except KeyError:
                    pass

    def cat_list(self, category: str) -> Set[T]:
        """Return all (active) objects from the specified 'categories'"""
        try:
            return self.categories[category]
        except KeyError:
            return set()

    def _create(self, name: str) -> T:
        """Create the object <T> from its name.
        Must be implemented in subclasses"""
        raise NotImplementedError

    def _initialize(self, obj: T) -> None:
        """Initialize the object <T> from its name.
        Must be implemented in subclasses"""
        raise NotImplementedError

    def _categorize(self, obj: T) -> List[str]:
        """List the categories of the object.
        Must be implemented in subclasses"""
        raise NotImplementedError
