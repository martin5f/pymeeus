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


from math import sqrt, degrees

from pyplanets.core.angle import Angle
from pyplanets.core.interpolation import Interpolation


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's now work with the Interpolation class
    print("\n" + 35 * "*")
    print("*** Use of Interpolation class")
    print(35 * "*" + "\n")

    i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])
    print("i = Interpolation([5, 3, 6, 1, 2, 4, 9], [10, 6, 12, 2, 4, 8])")
    print(i)
    print("NOTE:")
    print("   a. They are ordered in 'x'")
    print("   b. The extra value in 'x' was dropped")

    print("")

    # Use the copy constructor
    print("We can easily make a copy of an Interpolation object")
    j = Interpolation(i)
    print("j = Interpolation(i)")
    print(j)

    print("")

    # Get the number of interpolation points stored
    print_me("Number or interpolation points in 'j'", len(j))  # 6

    print("")

    j = Interpolation([0.0, 1.0, 3.0], [-1.0, -2.0, 2.0])
    print("j = Interpolation([0.0, 1.0, 3.0], [-1.0, -2.0, 2.0])")
    print(j)
    print_me("j(2)", j(2))
    print_me("j(0.5)", j(0.5))
    # Test with a value already in the data table
    print_me("j(1)", j(1))

    print("")

    # We can interpolate Angles too
    k = Interpolation(
        [27.0, 27.5, 28.0, 28.5, 29.0],
        [
            Angle(0, 54, 36.125),
            Angle(0, 54, 24.606),
            Angle(0, 54, 15.486),
            Angle(0, 54, 8.694),
            Angle(0, 54, 4.133),
        ],
    )
    print(
        "k = Interpolation([27.0, 27.5, 28.0, 28.5, 29.0],\n\
                      [Angle(0, 54, 36.125), Angle(0, 54, 24.606),\n\
                       Angle(0, 54, 15.486), Angle(0, 54, 8.694),\n\
                       Angle(0, 54, 4.133)])"
    )

    print_me("k(28.278)", Angle(k(28.278)).dms_str())

    print("")

    m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])
    print("m = Interpolation([-1.0, 0.0, 1.0], [-2.0, 3.0, 2.0])")
    print(m)
    # Get interpolated values
    print_me("m(-0.5)", m(-0.5))
    print_me("m(0.5)", m(0.5))
    # Get derivatives
    print_me("m'(-1.0)", m.derivative(-1.0))
    print_me("m'(-0.5)", m.derivative(-0.5))
    print_me("m'(0.0)", m.derivative(0.0))
    print_me("m'(0.5)", m.derivative(0.5))
    print_me("m'(1.0)", m.derivative(1.0))
    # Get the root within the interval
    print_me("m.root()", m.root())
    # Get the extremum within the interval
    print_me("m.minmax()", m.minmax())

    m = Interpolation(
        [29.43, 30.97, 27.69, 28.11, 31.58, 33.05],
        [
            0.4913598528,
            0.5145891926,
            0.4646875083,
            0.4711658342,
            0.5236885653,
            0.5453707057,
        ],
    )

    print_me("sin(29.5)\t", m(29.5))
    print_me("sin(30.0)\t", m(30.0))
    print_me("sin(30.5)\t", m(30.5))
    # Derivative must be adjusted because degrees were used instead of radians
    print_me("sin'(29.5)\t", degrees(m.derivative(29.5)))
    print_me("sin'(30.0)\t", degrees(m.derivative(30.0)))
    print_me("sqrt(3.0)/2.0\t", sqrt(3.0) / 2.0)
    print_me("sin'(30.5)\t", degrees(m.derivative(30.5)))


if __name__ == "__main__":
    main()
