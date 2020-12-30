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
from pyplanets.sun import Sun


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Sun functions
    print("\n" + 35 * "*")
    print("*** Use of Sun class")
    print(35 * "*" + "\n")

    # Compute an approximation of the Sun's true longitude
    epoch = Epoch(1992, 10, 13)
    true_lon, r = Sun.true_longitude_coarse(epoch)
    print_me("Sun's approximate true longitude", true_lon.dms_str(n_dec=0))
    # 199d 54' 36.0''
    print_me("Sun's radius vector", round(r, 5))  # 0.99766

    print("")

    # Now let's compute the Sun's approximate apparent longitude
    app_lon, r = Sun.apparent_longitude_coarse(epoch)
    print_me("Sun's approximate apparent longitude", app_lon.dms_str(n_dec=0))
    # 199d 54' 32.0''

    print("")

    # And now is the turn for the apparent right ascension and declination
    ra, delta, r = Sun.apparent_rightascension_declination_coarse(epoch)
    print_me("Sun's apparent right ascension", ra.ra_str(n_dec=1))
    # 13h 13' 31.4''
    print_me("Sun's apparent declination", delta.dms_str(n_dec=0))
    # -7d 47' 6.0''

    print("")

    # Let's compute Sun's postion, but more accurately
    epoch = Epoch(1992, 10, 13.0)
    l, b, r = Sun.geometric_geocentric_position(epoch, tofk5=False)
    print_me("Geometric Geocentric Longitude", round(l.to_positive(), 6))
    # 199.906016
    print_me("Geometric Geocentric Latitude", b.dms_str(n_dec=3))
    # 0.644''
    print_me("Radius vector", round(r, 8))
    # 0.99760775

    print("")

    # Compute Sun's apparent postion
    l, b, r = Sun.apparent_geocentric_position(epoch)
    print_me("Apparent Geocentric Longitude", l.to_positive().dms_str(n_dec=3))
    # 199d 54' 16.937''
    print_me("Apparent Geocentric Latitude", b.dms_str(n_dec=3))
    # 0.621''
    print_me("Radius vector", round(r, 8))
    # 0.99760775

    print("")

    # We can compute rectangular coordinates referred to mean equinox of date
    x, y, z = Sun.rectangular_coordinates_mean_equinox(epoch)
    print("Rectangular coordinates referred to mean equinox of date:")
    print_me("X", round(x, 7))  # -0.9379963
    print_me("Y", round(y, 6))  # -0.311654
    print_me("Z", round(z, 7))  # -0.1351207

    print("")

    # Now, compute rectangular coordinates w.r.t. standard equinox J2000.0
    x, y, z = Sun.rectangular_coordinates_j2000(epoch)
    print("Rectangular coordinates w.r.t. standard equinox J2000.0:")
    print_me("X", round(x, 8))  # -0.93740485
    print_me("Y", round(y, 8))  # -0.3131474
    print_me("Z", round(z, 8))  # -0.12456646

    print("")

    # Compute rectangular coordinates w.r.t. mean equinox of B1950.0
    x, y, z = Sun.rectangular_coordinates_b1950(epoch)
    print("Rectangular coordinates w.r.t. mean equinox of B1950.0:")
    print_me("X", round(x, 8))  # -0.94149557
    print_me("Y", round(y, 8))  # -0.30259922
    print_me("Z", round(z, 8))  # -0.11578695

    print("")

    # Compute rectangular coordinates w.r.t. an arbitrary mean equinox
    e_equinox = Epoch(2467616.0)
    x, y, z = Sun.rectangular_coordinates_equinox(epoch, e_equinox)
    print("Rectangular coordinates w.r.t. an arbitrary mean equinox:")
    print_me("X", round(x, 8))  # -0.93373777
    print_me("Y", round(y, 8))  # -0.32235109
    print_me("Z", round(z, 8))  # -0.12856709

    print("")

    # We can compute the date of equinoxes and solstices
    epoch = Sun.get_equinox_solstice(1962, target="summer")
    y, m, d, h, mi, s = epoch.get_full_date()
    print("The summer solstice of 1962:")
    print("{}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))
    # 1962/6/21 21:24:42.0

    epoch = Sun.get_equinox_solstice(2018, target="autumn")
    y, m, d, h, mi, s = epoch.get_full_date()
    print("The autumn equinox of 2018:")
    print("{}/{}/{} {}:{}:{}".format(y, m, d, h, mi, round(s, 0)))
    # 2018/9/23 1:55:14.0

    print("")

    # The equation of time, i.e., the difference between apparent and mean
    # time, can be easily computed
    epoch = Epoch(1992, 10, 13.0)
    m, s = Sun.equation_of_time(epoch)
    print("Equation of time difference: {} min {} secs".format(m, round(s, 1)))
    # 13m 42.6s

    print("")

    # Compute the ephemeris of physical observations of the Sun using
    # Carrington's formulas
    epoch = Epoch(1992, 10, 13)
    p, b0, l0 = Sun.ephemeris_physical_observations(epoch)
    print("Ephemeris of physical observations of the Sun:")
    print_me("P ", round(p, 2))  # 26.27
    print_me("B0", round(b0, 2))  # 5.99
    print_me("L0", round(l0, 2))  # 238.63

    print("")

    # Get the epoch when the Carrington's synodic rotation No. 'number' starts
    epoch = Sun.beginning_synodic_rotation(1699)
    print_me("Epoch for Carrington's synodic rotation No. 1699",
             round(epoch(), 3))  # 2444480.723


if __name__ == "__main__":
    main()
