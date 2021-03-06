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
from pyplanets.planets.mercury import Mercury


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Mercury class
    print("\n" + 35 * "*")
    print("*** Use of Mercury class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Mercury(epoch).geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    mercury = Mercury(epoch)
    earth = Earth(epoch)
    constellation = Constellation(earth, mercury)
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Mercury at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mercury(epoch).orbital_elements_mean_equinox()
    print_me("Mean longitude of the planet", round(l, 6))       # 203.494701
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))   # 0.38709831
    print_me("Eccentricity of the orbit", round(e, 7))          # 0.2056451
    print_me("Inclination on plane of the ecliptic", round(i, 6))   # 7.006171
    print_me("Longitude of the ascending node", round(ome, 5))  # 49.10765
    print_me("Argument of the perihelion", round(arg, 6))       # 29.367732

    print("")

    # Compute the time of the inferior conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury(epoch).inferior_conjunction()
    y, m, d = conjunction.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Inferior conjunction date", date)

    # Compute the time of the superior conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury(epoch).superior_conjunction()
    y, m, d = conjunction.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Superior conjunction date", date)

    print("")

    # Compute the time and angle of the western elongation close to 1993/11/1
    epoch = Epoch(1993, 11, 1.0)
    time, elongation = Mercury(epoch).western_elongation()
    y, m, d = time.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Western elongation date", date)
    elong = round(elongation, 4)
    print_me("Maximum western elongation angle", elong)

    print("")

    # Compute the time and angle of the eastern elongation close to 1990/8/1
    epoch = Epoch(1990, 8, 1.0)
    time, elongation = Mercury(epoch).eastern_elongation()
    y, m, d = time.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Eastern elongation date", date)
    elong = round(elongation, 4)
    print_me("Maximum eastern elongation angle", elong)

    print("")

    # Compute the time of the station in longitude #1 close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    sta1 = Mercury(epoch).station_longitude_1()
    y, m, d = sta1.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #1", date)

    # Compute the time of the station in longitude #2 close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    sta2 = Mercury(epoch).station_longitude_2()
    y, m, d = sta2.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #2", date)

    print("")

    # Find the epoch of the Perihelion closer to 2000/01/01
    epoch = Epoch(2000, 1, 1.0)
    e = Mercury(epoch).perihelion()
    y, m, d, h, mi, s = e.get_full_date()
    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'
    print_me("The Perihelion closest to 2000/1/1 happened on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    time, r = Mercury(epoch).passage_nodes()
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2018/11/24.7
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 0.3143


if __name__ == "__main__":
    main()
