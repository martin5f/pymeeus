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
from pymeeus_oo.calculation.Coordinates import (
    passage_nodes_elliptic
)
from pymeeus_oo.calculation.Epoch import Epoch
from pymeeus_oo.calculation.Interpolation import Interpolation
from pymeeus_oo.parameters.Uranus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.parameters.Venus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.planets.Planet import Planet

"""
.. module:: Uranus
   :synopsis: Class to model Uranus planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Uranus(Planet):
    """
    Class Uranus models that planet.
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
        >>> conj = Uranus.conjunction(epoch)
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1994
        >>> print(m)
        1
        >>> print(round(d, 4))
        12.7365
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Uranus' conjunction
        a = 2451579.489
        b = 369.656035
        m0 = 31.5219
        m1 = 4.333093
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute a couple auxiliary angles
        ee = 207.83 + 8.51 * t
        ff = 108.84 + 419.96 * t
        # Convert to radians
        ee = Angle(ee).rad()
        ff = Angle(ff).rad()
        corr = (-0.0859 + t * 0.0003 +
                sin(m) * (-3.8179 + t * (-0.0148 + t * 0.00003)) +
                cos(m) * (5.1228 + t * (-0.0105 - t * 0.00002)) +
                sin(2.0 * m) * (-0.0803 + t * 0.0011) +
                cos(2.0 * m) * (-0.1905 - t * 0.0006) +
                sin(3.0 * m) * (0.0088 + t * 0.0001) +
                cos(ee) * (0.885) +
                cos(ff) * (0.2153))
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

        >>> epoch = Epoch(1780, 12, 1.0)
        >>> oppo = Uranus.opposition(epoch)
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        1780
        >>> print(m)
        12
        >>> print(round(d, 4))
        17.5998
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Uranus' opposition
        a = 2451764.317
        b = 369.656035
        m0 = 213.6884
        m1 = 4.333093
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute a couple auxiliary angles
        ee = 207.83 + 8.51 * t
        ff = 108.84 + 419.96 * t
        # Convert to radians
        ee = Angle(ee).rad()
        ff = Angle(ff).rad()
        corr = (0.0844 - t * 0.0006 +
                sin(m) * (-0.1048 + t * 0.0246) +
                cos(m) * (-5.1221 + t * (0.0104 + t * 0.00003)) +
                sin(2.0 * m) * (-0.1428 + t * 0.0005) +
                cos(2.0 * m) * (-0.0148 - t * 0.0013) +
                cos(3.0 * m) * (0.0055) +
                cos(ee) * (0.885) +
                cos(ff) * (0.2153))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def perihelion_aphelion(epoch: Epoch, perihelion=True) -> Epoch:
        """This method computes the time of Perihelion (or Aphelion) closer to
        a given epoch.

        :param epoch: Epoch close to the desired Perihelion (or Aphelion)
        :type epoch: :py:class:`Epoch`
        :param peihelion: If True, the epoch of the closest Perihelion is
            computed, if False, the epoch of the closest Aphelion is found.
        :type bool:

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input values are of wrong type.

        .. note:: The solution provided by this method may have several days of
            error.

        >>> epoch = Epoch(1880, 1, 1.0)
        >>> e = Uranus.perihelion_aphelion(epoch)
        >>> y, m, d = e.get_date()
        >>> print(y)
        1882
        >>> print(m)
        3
        >>> print(int(d))
        18
        >>> epoch = Epoch(2090, 1, 1.0)
        >>> e = Uranus.perihelion_aphelion(epoch, perihelion=False)
        >>> y, m, d = e.get_date()
        >>> print(y)
        2092
        >>> print(m)
        11
        >>> print(int(d))
        22
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input value")
        # First approximation
        k = 0.0119 * (epoch.year() - 2051.1)
        if perihelion:
            k = round(k)
        else:
            k = round(k + 0.5) - 0.5
        jde = 2470213.5 + k * (30694.8767 - k * 0.00541)
        # Compute the epochs 1 year before and after
        jde_before = jde - 360.0
        jde_after = jde + 360.0
        # Compute the Sun-Uranus distance for each epoch
        l, b, r_b = Uranus.geometric_heliocentric_position(Epoch(jde_before))
        l, b, r = Uranus.geometric_heliocentric_position(Epoch(jde))
        l, b, r_a = Uranus.geometric_heliocentric_position(Epoch(jde_after))
        # Call an interpolation object
        m = Interpolation([jde_before, jde, jde_after], [r_b, r, r_a])
        sol = m.minmax()
        return Epoch(sol)

    @staticmethod
    def passage_nodes(epoch: Epoch, ascending=True) -> (Epoch, float):
        """This function computes the time of passage by the nodes (ascending
        or descending) of Uranus, nearest to the given epoch.

        :param epoch: Epoch closest to the node passage
        :type epoch: :py:class:`Epoch`
        :param ascending: Whether the time of passage by the ascending (True)
            or descending (False) node will be computed
        :type ascending: bool

        :returns: Tuple containing:
            - Time of passage through the node (:py:class:`Epoch`)
            - Radius vector when passing through the node (in AU, float)
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(2019, 1, 1)
        >>> time, r = Uranus.passage_nodes(epoch)
        >>> year, month, day = time.get_date()
        >>> print(year)
        2028
        >>> print(month)
        8
        >>> print(round(day, 1))
        23.2
        >>> print(round(r, 4))
        19.3201
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Get the orbital parameters
        l, a, e, i, ome, arg = Uranus.orbital_elements_mean_equinox(epoch)
        # Compute the time of passage through perihelion
        t = Uranus.perihelion_aphelion(epoch)
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r

    @staticmethod
    def magnitude(sun_dist, earth_dist):
        """This function computes the approximate magnitude of Uranus.

        :param sun_dist: Distance from Uranus to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Uranus to Earth, in Astronomical Units
        :type earth_dist: float

        :returns: Uranus's magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.
        """

        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)):
            raise TypeError("Invalid input types")
        m = -6.85 + 5.0 * log10(sun_dist * earth_dist)
        return round(m, 1)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Uranus class
    print("\n" + 35 * "*")
    print("*** Use of Uranus class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Uranus.geometric_heliocentric_position(epoch)
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    ra, dec, elon = Uranus.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Uranus at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Uranus.orbital_elements_mean_equinox(epoch)
    print_me("Mean longitude of the planet", round(l, 6))  # 235.517526
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))  # 19.21844604
    print_me("Eccentricity of the orbit", round(e, 7))  # 0.0463634
    print_me("Inclination on plane of the ecliptic", round(i, 6))  # 0.77372
    print_me("Longitude of the ascending node", round(ome, 5))  # 74.34776
    print_me("Argument of the perihelion", round(arg, 6))  # 99.630865

    print("")

    # Compute the time of the conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conj = Uranus.conjunction(epoch)
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 1780/12/1
    epoch = Epoch(1780, 12, 1.0)
    oppo = Uranus.opposition(epoch)
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)

    print("")

    # Find the epoch of the Perihelion closer to 1780/1/1
    epoch = Epoch(1780, 1, 1.0)
    e = Uranus.perihelion_aphelion(epoch)
    y, m, d = e.get_date()
    peri = str(y) + '/' + str(m) + '/' + str(int(d))
    print_me("The Perihelion closest to 1780/1/1 happened on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    time, r = Uranus.passage_nodes(epoch)
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2028/8/23.2
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 19.3201


if __name__ == "__main__":

    main()
