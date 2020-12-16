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
from pymeeus_oo.parameters.Mercury_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.parameters.Venus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.planets.Planet import Planet

"""
.. module:: Mercury
   :synopsis: Class to model Mercury planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Mercury(Planet):
    """
    Class Mercury models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    @staticmethod
    def inferior_conjunction(epoch: Epoch) -> Epoch:
        """This method computes the time of the inferior conjunction closest to
        the given epoch.

        :param epoch: Epoch close to the desired inferior conjunction
        :type epoch: :py:class:`Epoch`

        :returns: The time when the inferior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conjunction = Mercury.inferior_conjunction(epoch)
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1993
        >>> print(m)
        11
        >>> print(round(d, 4))
        6.1449
        >>> epoch = Epoch(1631, 10, 1.0)
        >>> conjunction = Mercury.inferior_conjunction(epoch)
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1631
        >>> print(m)
        11
        >>> print(round(d, 3))
        7.306
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's inferior conjunction
        a = 2451612.023
        b = 115.8774771
        m0 = 63.5867
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (0.0545 + 0.0002 * t +
                sin(m) * (-6.2008 + t * (0.0074 + t * 0.00003)) +
                cos(m) * (-3.275 + t * (-0.0197 + t * 0.00001)) +
                sin(2.0 * m) * (0.4737 + t * (-0.0052 - t * 0.00001)) +
                cos(2.0 * m) * (0.8111 + t * (0.0033 - t * 0.00002)) +
                sin(3.0 * m) * (0.0037 + t * 0.0018) +
                cos(3.0 * m) * (-0.1768 + t * t * 0.00001) +
                sin(4.0 * m) * (-0.0211 - t * 0.0004) +
                cos(4.0 * m) * (0.0326 - t * 0.0003) +
                sin(5.0 * m) * (0.0083 + t * 0.0001) +
                cos(5.0 * m) * (-0.004 + t * 0.0001))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def superior_conjunction(epoch: Epoch) -> Epoch:
        """This method computes the time of the superior conjunction closest to
        the given epoch.

        :param epoch: Epoch close to the desired superior conjunction
        :type epoch: :py:class:`Epoch`

        :returns: The time when the superior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conjunction = Mercury.superior_conjunction(epoch)
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1993
        >>> print(m)
        8
        >>> print(round(d, 4))
        29.3301
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's superior conjunction
        a = 2451554.084
        b = 115.8774771
        m0 = 6.4822
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-0.0548 - 0.0002 * t +
                sin(m) * (7.3894 + t * (-0.01 - t * 0.00003)) +
                cos(m) * (3.22 + t * (0.0197 - t * 0.00001)) +
                sin(2.0 * m) * (0.8383 + t * (-0.0064 - t * 0.00001)) +
                cos(2.0 * m) * (0.9666 + t * (0.0039 - t * 0.00003)) +
                sin(3.0 * m) * (0.077 - t * 0.0026) +
                cos(3.0 * m) * (0.2758 + t * (0.0002 - t * 0.00002)) +
                sin(4.0 * m) * (-0.0128 - t * 0.0008) +
                cos(4.0 * m) * (0.0734 + t * (-0.0004 - t * 0.00001)) +
                sin(5.0 * m) * (-0.0122 - t * 0.0002) +
                cos(5.0 * m) * (0.0173 - t * 0.0002))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def western_elongation(epoch: Epoch) -> Epoch:
        """This method computes the time of the western elongation closest to
        the given epoch, as well as the corresponding maximum elongation angle.

        :param epoch: Epoch close to the desired western elongation
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the time when the western elongation happens, as
            an Epoch, and the maximum elongation angle, as an Angle
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 11, 1.0)
        >>> time, elongation = Mercury.western_elongation(epoch)
        >>> y, m, d = time.get_date()
        >>> print(y)
        1993
        >>> print(m)
        11
        >>> print(round(d, 4))
        22.6386
        >>> print(round(elongation, 4))
        19.7506
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's inferior conjunction
        a = 2451612.023
        b = 115.8774771
        m0 = 63.5867
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (21.6249 - 0.0002 * t +
                sin(m) * (0.1306 + t * 0.0065) +
                cos(m) * (-2.7661 + t * (-0.0011 + t * 0.00001)) +
                sin(2.0 * m) * (0.2438 + t * (-0.0024 - t * 0.00001)) +
                cos(2.0 * m) * (0.5767 + t * 0.0023) +
                sin(3.0 * m) * (0.1041) +
                cos(3.0 * m) * (-0.0184 + t * 0.0007) +
                sin(4.0 * m) * (-0.0051 - t * 0.0001) +
                cos(4.0 * m) * (0.0048 + t * 0.0001) +
                sin(5.0 * m) * (0.0026) +
                cos(5.0 * m) * (0.0037))
        elon = (22.4143 - 0.0001 * t +
                sin(m) * (4.3651 + t * (-0.0048 - t * 0.00002)) +
                cos(m) * (2.3787 + t * (0.0121 - t * 0.00001)) +
                sin(2.0 * m) * (0.2674 + t * 0.0022) +
                cos(2.0 * m) * (-0.3873 + t * (0.0008 + t * 0.00001)) +
                sin(3.0 * m) * (-0.0369 - t * 0.0001) +
                cos(3.0 * m) * (0.0017 - t * 0.0001) +
                sin(4.0 * m) * (0.0059) +
                cos(4.0 * m) * (0.0061 + t * 0.0001) +
                sin(5.0 * m) * (0.0007) +
                cos(5.0 * m) * (-0.0011))
        elon = Angle(elon).to_positive()
        to_return = jde0 + corr
        return Epoch(to_return), elon

    @staticmethod
    def eastern_elongation(epoch: Epoch) -> Epoch:
        """This method computes the time of the eastern elongation closest to
        the given epoch, as well as the corresponding maximum elongation angle.

        :param epoch: Epoch close to the desired eastern elongation
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the time when the eastern elongation happens, as
            an Epoch, and the maximum elongation angle, as an Angle
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1990, 8, 1.0)
        >>> time, elongation = Mercury.eastern_elongation(epoch)
        >>> y, m, d = time.get_date()
        >>> print(y)
        1990
        >>> print(m)
        8
        >>> print(round(d, 4))
        11.8514
        >>> print(round(elongation, 4))
        27.4201
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's inferior conjunction
        a = 2451612.023
        b = 115.8774771
        m0 = 63.5867
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-21.6101 + 0.0002 * t +
                sin(m) * (-1.9803 + t * (-0.006 + t * 0.00001)) +
                cos(m) * (1.4151 + t * (-0.0072 - t * 0.00001)) +
                sin(2.0 * m) * (0.5528 + t * (-0.0005 - t * 0.00001)) +
                cos(2.0 * m) * (0.2905 + t * (0.0034 + t * 0.00001)) +
                sin(3.0 * m) * (-0.1121 + t * (-0.0001 + t * 0.00001)) +
                cos(3.0 * m) * (-0.0098 - t * 0.0015) +
                sin(4.0 * m) * (0.0192) +
                cos(4.0 * m) * (0.0111 + t * 0.0004) +
                sin(5.0 * m) * (-0.0061) +
                cos(5.0 * m) * (-0.0032 - t * 0.0001))
        elon = (22.4697 +
                sin(m) * (-4.2666 + t * (0.0054 + t * 0.00002)) +
                cos(m) * (-1.8537 - t * 0.0137) +
                sin(2.0 * m) * (0.3598 + t * (0.0008 - t * 0.00001)) +
                cos(2.0 * m) * (-0.068 + t * 0.0026) +
                sin(3.0 * m) * (-0.0524 - t * 0.0003) +
                cos(3.0 * m) * (0.0052 - t * 0.0006) +
                sin(4.0 * m) * (0.0107 + t * 0.0001) +
                cos(4.0 * m) * (-0.0013 + t * 0.0001) +
                sin(5.0 * m) * (-0.0021) +
                cos(5.0 * m) * (0.0003))
        elon = Angle(elon).to_positive()
        to_return = jde0 + corr
        return Epoch(to_return), elon

    @staticmethod
    def station_longitude_1(epoch: Epoch) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired inferior conjunction
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 1st statin in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> sta1 = Mercury.station_longitude_1(epoch)
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        1993
        >>> print(m)
        10
        >>> print(round(d, 4))
        25.9358
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's inferior conjunction
        a = 2451612.023
        b = 115.8774771
        m0 = 63.5867
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-11.0761 + 0.0003 * t +
                sin(m) * (-4.7321 + t * (0.0023 + t * 0.00002)) +
                cos(m) * (-1.323 - t * 0.0156) +
                sin(2.0 * m) * (0.227 - t * 0.0046) +
                cos(2.0 * m) * (0.7184 + t * (0.0013 - t * 0.00002)) +
                sin(3.0 * m) * (0.0638 + t * 0.0016) +
                cos(3.0 * m) * (-0.1655 + t * 0.0007) +
                sin(4.0 * m) * (-0.0395 - t * 0.0003) +
                cos(4.0 * m) * (0.0247 - t * 0.0006) +
                sin(5.0 * m) * (0.0131) +
                cos(5.0 * m) * (-0.0008 + t * 0.0002))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def station_longitude_2(epoch: Epoch) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired inferior conjunction
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> sta2 = Mercury.station_longitude_2(epoch)
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        1993
        >>> print(m)
        11
        >>> print(round(d, 4))
        15.0724
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mercury's inferior conjunction
        a = 2451612.023
        b = 115.8774771
        m0 = 63.5867
        m1 = 114.2088742
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (11.1343 - 0.0001 * t +
                sin(m) * (-3.9137 + t * (0.0073 + t * 0.00002)) +
                cos(m) * (-3.3861 + t * (-0.0128 + t * 0.00001)) +
                sin(2.0 * m) * (0.5222 + t * (-0.004 - t * 0.00002)) +
                cos(2.0 * m) * (0.5929 + t * (0.0039 - t * 0.00002)) +
                sin(3.0 * m) * (-0.0593 + t * 0.0018) +
                cos(3.0 * m) * (-0.1733 * t * (-0.0007 + t * 0.00001)) +
                sin(4.0 * m) * (-0.0053 - t * 0.0006) +
                cos(4.0 * m) * (0.0476 - t * 0.0001) +
                sin(5.0 * m) * (0.007 + t * 0.0002) +
                cos(5.0 * m) * (-0.0115 + t * 0.0001))
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

        >>> epoch = Epoch(2000, 1, 1.0)
        >>> e = Mercury.perihelion_aphelion(epoch)
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2000
        >>> print(m)
        2
        >>> print(d)
        15
        >>> print(h)
        18
        >>> epoch = Epoch(2000, 3, 1.0)
        >>> e = Mercury.perihelion_aphelion(epoch, perihelion=False)
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2000
        >>> print(m)
        3
        >>> print(d)
        30
        >>> print(h)
        17
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input value")
        # First approximation
        k = 4.15201 * (epoch.year() - 2000.12)
        if perihelion:
            k = round(k)
        else:
            k = round(k + 0.5) - 0.5
        jde = 2451590.257 + k * 87.96934963
        # Compute the epochs half a day before and after
        jde_before = jde - 0.5
        jde_after = jde + 0.5
        # Compute the Sun-Mercury distance for each epoch
        l, b, r_b = Mercury.geometric_heliocentric_position(Epoch(jde_before))
        l, b, r = Mercury.geometric_heliocentric_position(Epoch(jde))
        l, b, r_a = Mercury.geometric_heliocentric_position(Epoch(jde_after))
        # Call an interpolation object
        m = Interpolation([jde_before, jde, jde_after], [r_b, r, r_a])
        sol = m.minmax()
        return Epoch(sol)

    @staticmethod
    def passage_nodes(epoch: Epoch, ascending=True) -> (Epoch, float):
        """This function computes the time of passage by the nodes (ascending
        or descending) of Mercury, nearest to the given epoch.

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
        >>> time, r = Mercury.passage_nodes(epoch)
        >>> year, month, day = time.get_date()
        >>> print(year)
        2018
        >>> print(month)
        11
        >>> print(round(day, 1))
        24.7
        >>> print(round(r, 4))
        0.3143
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Get the orbital parameters
        l, a, e, i, ome, arg = Mercury.orbital_elements_mean_equinox(epoch)
        # Compute the time of passage through perihelion
        t = Mercury.perihelion_aphelion(epoch)
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r

    @staticmethod
    def magnitude(sun_dist, earth_dist, phase_angle):
        """This function computes the approximate magnitude of Mercury.

        :param sun_dist: Distance from Mercury to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance Mercury to Earth, in Astronomical Units
        :type earth_dist: float
        :param phase_angle: Mercury phase angle
        :type phase_angle: float, :py:class:`Angle`

        :returns: Mercury's magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.
        """

        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)
                and isinstance(phase_angle, (float, Angle))):
            raise TypeError("Invalid input types")
        i = float(phase_angle)
        i50 = i - 50.0
        m = (1.16 + 5.0 * log10(sun_dist * earth_dist) + 0.02838 * i50 +
             0.0001023 * i50 * i50)
        return round(m, 1)


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
    lon, lat, r = Mercury.geometric_heliocentric_position(epoch)
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    ra, dec, elon = Mercury.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Mercury at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Mercury.orbital_elements_mean_equinox(epoch)
    print_me("Mean longitude of the planet", round(l, 6))       # 203.494701
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))   # 0.38709831
    print_me("Eccentricity of the orbit", round(e, 7))          # 0.2056451
    print_me("Inclination on plane of the ecliptic", round(i, 6))   # 7.006171
    print_me("Longitude of the ascending node", round(ome, 5))  # 49.10765
    print_me("Argument of the perihelion", round(arg, 6))       # 29.367732

    print("")

    # Compute the time of the inferior conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury.inferior_conjunction(epoch)
    y, m, d = conjunction.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Inferior conjunction date", date)

    # Compute the time of the superior conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conjunction = Mercury.superior_conjunction(epoch)
    y, m, d = conjunction.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Superior conjunction date", date)

    print("")

    # Compute the time and angle of the western elongation close to 1993/11/1
    epoch = Epoch(1993, 11, 1.0)
    time, elongation = Mercury.western_elongation(epoch)
    y, m, d = time.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Western elongation date", date)
    elong = round(elongation, 4)
    print_me("Maximum western elongation angle", elong)

    print("")

    # Compute the time and angle of the eastern elongation close to 1990/8/1
    epoch = Epoch(1990, 8, 1.0)
    time, elongation = Mercury.eastern_elongation(epoch)
    y, m, d = time.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Eastern elongation date", date)
    elong = round(elongation, 4)
    print_me("Maximum eastern elongation angle", elong)

    print("")

    # Compute the time of the station in longitude #1 close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    sta1 = Mercury.station_longitude_1(epoch)
    y, m, d = sta1.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #1", date)

    # Compute the time of the station in longitude #2 close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    sta2 = Mercury.station_longitude_2(epoch)
    y, m, d = sta2.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #2", date)

    print("")

    # Find the epoch of the Perihelion closer to 2000/01/01
    epoch = Epoch(2000, 1, 1.0)
    e = Mercury.perihelion_aphelion(epoch)
    y, m, d, h, mi, s = e.get_full_date()
    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'
    print_me("The Perihelion closest to 2000/1/1 happened on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    time, r = Mercury.passage_nodes(epoch)
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2018/11/24.7
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 0.3143


if __name__ == "__main__":

    main()
