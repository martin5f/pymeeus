# -*- coding: utf-8 -*-


# PyMeeus: Python module implementing astronomical algorithms.
# Copyright (C) 2018  Dagoberto Salazar
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


from math import sin, cos, log10

from pymeeus_oo.calculation.Angle import Angle
from pymeeus_oo.calculation.Epoch import Epoch
from pymeeus_oo.parameters.Neptune_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.planets.Planet import Planet

"""
.. module:: Neptune
   :synopsis: Class to model Neptune planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Neptune(Planet):
    """
    Class Neptune models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    @staticmethod
    def conjunction(epoch: Epoch) -> Epoch:
        """This method computes the time of the conjunction closest to the
        given epoch.

        :param epoch: Epoch close to the desired conjunction
        :type epoch: :py:class:`Epoch`

        :returns: The time when the conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conj = Neptune.conjunction(epoch)
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1994
        >>> print(m)
        1
        >>> print(round(d, 4))
        11.3057
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Neptune's conjunction
        a = 2451569.379
        b = 367.486703
        m0 = 21.5569
        m1 = 2.194998
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute a couple auxiliary angles
        ee = 207.83 + 8.51 * t
        gg = 276.74 + 209.98 * t
        # Convert to radians
        ee = Angle(ee).rad()
        gg = Angle(gg).rad()
        corr = (0.0168 +
                sin(m) * (-2.5606 + t * (0.0088 + t * 0.00002)) +
                cos(m) * (-0.8611 + t * (-0.0037 + t * 0.00002)) +
                sin(2.0 * m) * (0.0118 + t * (-0.0004 + t * 0.00001)) +
                cos(2.0 * m) * (0.0307 - t * 0.0003) +
                cos(ee) * (-0.5964) +
                cos(gg) * (0.0728))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def opposition(epoch: Epoch) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1846, 8, 1)
        >>> oppo = Neptune.opposition(epoch)
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        1846
        >>> print(m)
        8
        >>> print(round(d, 4))
        20.1623
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Neptune's opposition
        a = 2451753.122
        b = 367.486703
        m0 = 202.6544
        m1 = 2.194998
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute a couple auxiliary angles
        ee = 207.83 + 8.51 * t
        gg = 276.74 + 209.98 * t
        # Convert to radians
        ee = Angle(ee).rad()
        gg = Angle(gg).rad()
        corr = (-0.014 + t * t * 0.00001 +
                sin(m) * (-1.3486 + t * (0.001 + t * 0.00001)) +
                cos(m) * (0.8597 + t * 0.0037) +
                sin(2.0 * m) * (-0.0082 + t * (-0.0002 + t * 0.00001)) +
                cos(2.0 * m) * (0.0037 - t * 0.0003) +
                cos(ee) * (-0.5964) +
                cos(gg) * (0.0728))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def magnitude(sun_dist, earth_dist):
        """This function computes the approximate magnitude of Neptune.

        :param sun_dist: Distance from Neptune to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance Neptune to Earth, in Astronomical Units
        :type earth_dist: float

        :returns: Neptune's magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.
        """

        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)):
            raise TypeError("Invalid input types")
        m = -7.05 + 5.0 * log10(sun_dist * earth_dist)
        return round(m, 1)


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
    lon, lat, r = Neptune.geometric_heliocentric_position(epoch)
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    ra, dec, elon = Neptune.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Neptune at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Neptune.orbital_elements_mean_equinox(epoch)
    print_me("Mean longitude of the planet", round(l, 6))  # 88.321947
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))  # 30.11038676
    print_me("Eccentricity of the orbit", round(e, 7))  # 0.0094597
    print_me("Inclination on plane of the ecliptic", round(i, 6))  # 1.763855
    print_me("Longitude of the ascending node", round(ome, 5))  # 132.46986
    print_me("Argument of the perihelion", round(arg, 6))  # -83.415521

    print("")

    # Compute the time of the conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conj = Neptune.conjunction(epoch)
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 1846/8/1
    epoch = Epoch(1846, 8, 1)
    oppo = Neptune.opposition(epoch)
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)


if __name__ == "__main__":

    main()
