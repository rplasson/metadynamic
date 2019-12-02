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


from typing import Any


class InvalidError(ValueError):
    pass


class Invalid:
    _invalrepr = "Invalid Object"
    _invalstr = ""

    def __bool__(self) -> bool:
        return False

    def __str__(self) -> str:
        return self._invalstr

    def __repr__(self) -> str:
        return self._invalrepr


def isvalid(obj: Any) -> bool:
    return not isinstance(obj, Invalid)


class InvalidInt(Invalid, int):
    _invalrepr = "Invalid int"


class InvalidFloat(Invalid, float):
    _invalrepr = "Invalid float"


class InvalidStr(Invalid, str):
    _invalrepr = "Invalid string"


def invalfactory(baseclass, name):
    class InvalidClass(Invalid, baseclass):
        _invalrepr = name

    return InvalidClass

