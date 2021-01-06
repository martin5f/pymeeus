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


from examples.pymeeus import usage_venus, usage_angle, usage_neptune, usage_minor, usage_mars, usage_pluto, \
    usage_uranus, usage_saturn, usage_epoch, usage_base, usage_earth, usage_interpolation, usage_curvefitting, \
    usage_coordinates, usage_sun, usage_jupiter, usage_mercury


def main():
    # Base
    usage_angle.main()
    usage_base.main()
    usage_coordinates.main()
    usage_curvefitting.main()
    usage_epoch.main()
    usage_interpolation.main()

    # Sun
    usage_sun.main()

    # Planets
    usage_mercury.main()
    usage_venus.main()
    usage_earth.main()
    usage_mars.main()
    usage_jupiter.main()
    usage_saturn.main()
    usage_uranus.main()
    usage_neptune.main()
    usage_pluto.main()
    usage_minor.main()


if __name__ == "__main__":
    main()
