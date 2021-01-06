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

from pyplanets.planets.neptune import Neptune


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Neptune class
    print("\n" + 35 * "*")
    print("*** Use of Neptune class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Neptune(epoch).geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    constellation = Constellation(Earth(epoch), Neptune(epoch))
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Neptune at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Neptune(epoch).orbital_elements_mean_equinox()
    print_me("Mean longitude of the planet", round(l, 6))  # 88.321947
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))  # 30.11038676
    print_me("Eccentricity of the orbit", round(e, 7))  # 0.0094597
    print_me("Inclination on plane of the ecliptic", round(i, 6))  # 1.763855
    print_me("Longitude of the ascending node", round(ome, 5))  # 132.46986
    print_me("Argument of the perihelion", round(arg, 6))  # -83.415521

    print("")

    # Compute the time of the conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conj = Neptune(epoch).conjunction()
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 1846/8/1
    epoch = Epoch(1846, 8, 1)
    oppo = Neptune(epoch).opposition()
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)


if __name__ == "__main__":
    main()
