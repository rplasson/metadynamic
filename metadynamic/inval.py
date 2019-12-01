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


from typing import Any, Generic, TypeVar

T = TypeVar("T")


class Invalid(Generic[T]):
    _invalrepr = "Invalid Object"

    def __bool__(self) -> bool:
        return False

    def __str__(self) -> str:
        return self._invalrepr

    def __repr__(self) -> str:
        return self._invalrepr

    def validate(self) -> T:
        raise NotImplementedError


def isvalid(obj: Any) -> bool:
    if isinstance(obj, Invalid):
        return False
    return True


class InvalInt(Invalid[int], int):
    _invalrepr = "Invalid int"
    def validate(self) -> int:
        return int(self)


class InvalFloat(Invalid[float], float):
    _invalrepr = "Invalid float"


invalint = InvalInt()
invalfloat = InvalFloat()


def testint(a: int = invalint) -> None:
    assert isinstance(a, int)
    if isvalid(a):
        print(a)
    else:
        print("Give a value")


def testfloat(a: float = invalfloat) -> None:
    assert isinstance(a, float)
    if isvalid(a):
        print(a)
    else:
        print("Give a value")
