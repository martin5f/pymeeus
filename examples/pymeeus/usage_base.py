# -*- coding: utf-8 -*-


# PyPlanets: Object-oriented refactoring of PyMeeus, a Python library implementing astronomical algorithms.
# Copyright (C) 2020  Martin Fünffinger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


from pyplanets.core.base import TOL, machine_accuracy
from pyplanets.core.base import get_ordinal_suffix


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's print the tolerance
    print_me("The default value for the tolerance is", TOL)

    # Find the accuracy of this computer
    j, d = machine_accuracy()
    print_me("Number of significant BITS in the mantissa\t", j)
    print_me("Number of significant DIGITS in a decimal number", d)

    print("")

    print_me("The suffix for ordinal 2 is", get_ordinal_suffix(2))
    print_me("The suffix for ordinal 11 is", get_ordinal_suffix(11))
    print_me("The suffix for ordinal 12 is", get_ordinal_suffix(12))
    print_me("The suffix for ordinal 13 is", get_ordinal_suffix(13))
    print_me("The suffix for ordinal 14 is", get_ordinal_suffix(14))
    print_me("The suffix for ordinal 16 is", get_ordinal_suffix(16))
    print_me("The suffix for ordinal 23 is", get_ordinal_suffix(23))


if __name__ == "__main__":
    main()
