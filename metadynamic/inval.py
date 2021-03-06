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

"""Define invalid values of any type.

Invalid objects of type X can be obtained by deriving an "Invalid X" class from both X and
L{Invalid}. Objects validity can be used with the L{isvalidf} function. L{InvalidError} shall be
used when an invalid value is found where a valid one is expected.  L{invalidint}, L{invalidfloat}
and L{invalidstr} are defined as examples.

"""

from typing import Any


class InvalidError(ValueError):
    """Error to be raised when an invalid value is found where unexpected."""


class Invalid:
    """Define Invalid classes of any type.

    An invalid object of type <T> is an instance of a class
    that inherits from both Invalid and <T>

        >>> class InvalidInt(Invalid, int):
                _invalrepr = "Invalid int"
        >>> invalidint = InvalidInt()
        >>> (invalidint, 3)
        (Invalid int, 3)
        >>> isvalid(invalidint), isvalid(3)
        (False, True)
        >>> isinstance(invalidint, int), isinstance(3,int)
        (True, True)

    Can be useful to define a default "no value"
    parameter, that will be seen as an instance of a regular type

    _invalrepr and _invalstr shall be overriden in subclasses

    """

    _invalrepr = "Invalid Object"
    """Representation of the invalid object"""
    _invalstr = ""
    """String conversion of the invalid object"""

    def __bool__(self) -> bool:
        """Return False from any invalid object."""
        return False

    def __str__(self) -> str:
        """Convert invalid objects to self._invalstr.

        defaults to empty string, shall be overriden in subclasses if necessary.

        """
        return self._invalstr

    def __repr__(self) -> str:
        """Represent invalid objects by self._invalrepr.

        defaults to 'Invalid Object', shall be overriden in subclasses.

        """
        return self._invalrepr


def isvalid(obj: Any) -> bool:
    """Return False only if the object is an instance of Invalid class.

    @param obj: object to be tested
    @return: is the object valid?
    @rtype: bool

    """
    return not isinstance(obj, Invalid)


class InvalidInt(Invalid, int):
    """Invalid int class."""

    _invalrepr = "Invalid int"


invalidint: InvalidInt = InvalidInt()
"""invalint int object"""


class InvalidFloat(Invalid, float):
    """Invalit float class."""

    _invalrepr = "Invalid float"


invalidfloat: InvalidFloat = InvalidFloat()
"""invalid float object"""


class InvalidStr(Invalid, str):
    """invalid str class."""

    _invalrepr = "Invalid string"


invalidstr: InvalidStr = InvalidStr()
"""invalid str object"""
