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


from pyplanets.core.epoch import Epoch
from pyplanets.core.constellation import Constellation
from pyplanets.planets.earth import Earth

from pyplanets.planets.pluto import Pluto


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Pluto class
    print("\n" + 35 * "*")
    print("*** Use of Pluto class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = Pluto(epoch).geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/10/13:
    epoch = Epoch(1992, 10, 13.0)
    ra, dec = Pluto.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))

    constellation = Constellation(Earth(epoch), Pluto(epoch))
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

if __name__ == "__main__":
    main()
