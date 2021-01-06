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


from pymeeus import Earth, Angle, base, Coordinates, CurveFitting, Epoch, Interpolation, Sun
from pymeeus import Jupiter
from pymeeus import Mars
from pymeeus import Mercury
from pymeeus import Minor
from pymeeus import Neptune
from pymeeus import Pluto
from pymeeus import Saturn
from pymeeus import Uranus
from pymeeus import Venus


def main():
    # Base
    Angle.main()
    base.main()
    Coordinates.main()
    CurveFitting.main()
    Epoch.main()
    Interpolation.main()

    # Sun
    Sun.main()

    # Planets
    Mercury.main()
    Venus.main()
    Earth.main()
    Mars.main()
    Jupiter.main()
    Saturn.main()
    Uranus.main()
    Neptune.main()
    Pluto.main()
    Minor.main()


if __name__ == "__main__":
    main()
