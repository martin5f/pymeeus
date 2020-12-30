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


from math import sin, cos, log10

from pyplanets.core.angle import Angle
from pyplanets.core.epoch import Epoch
from pyplanets.parameters.saturn_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Saturn
   :synopsis: Class to model Saturn planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Saturn(Planet):
    """
    Class Saturn models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    def conjunction(self) -> Epoch:
        """This method computes the time of the conjunction closest to the
        given epoch.

        :returns: The time when the conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2125, 6, 1.0)
        >>> conj = Saturn(epoch).conjunction()
        >>> y, m, d = conj.get_date()
        >>> print(y)
        2125
        >>> print(m)
        8
        >>> print(round(d, 4))
        26.4035
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's conjunction
        a = 2451681.124
        b = 378.091904
        m0 = 131.6934
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute auxiliary angles
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (0.0172 + t * (-0.0006 + t * 0.00023) +
                sin(m) * (-8.5885 + t * (0.0411 + t * 0.0002)) +
                cos(m) * (-1.147 + t * (0.0352 - t * 0.00011)) +
                sin(2.0 * m) * (0.3331 + t * (-0.0034 - t * 0.00001)) +
                cos(2.0 * m) * (0.1145 + t * (-0.0045 + t * 0.00002)) +
                sin(3.0 * m) * (-0.0169 + t * 0.0002) +
                cos(3.0 * m) * (-0.0109 + t * 0.0004) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def opposition(self) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(-6, 9, 1.0)
        >>> oppo = Saturn(epoch).opposition()
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        -6
        >>> print(m)
        9
        >>> print(round(d, 4))
        14.3709
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (-0.0209 + t * (0.0006 + t * 0.00023) +
                sin(m) * (4.5795 + t * (-0.0312 - t * 0.00017)) +
                cos(m) * (1.1462 + t * (-0.0351 + t * 0.00011)) +
                sin(2.0 * m) * (0.0985 - t * 0.0015) +
                cos(2.0 * m) * (0.0733 + t * (-0.0031 + t * 0.00001)) +
                sin(3.0 * m) * (0.0025 - t * 0.0001) +
                cos(3.0 * m) * (0.005 - t * 0.0002) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_1(self) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 11, 1.0)
        >>> sta1 = Saturn(epoch).station_longitude_1()
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        2018
        >>> print(m)
        4
        >>> print(round(d, 4))
        17.9433
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (-68.884 + t * (0.0009 + t * 0.00023) +
                sin(m) * (5.5452 + t * (-0.0279 - t * 0.0002)) +
                cos(m) * (3.0727 + t * (-0.043 + t * 0.00007)) +
                sin(2.0 * m) * (0.1101 + t * (-0.0006 - t * 0.00001)) +
                cos(2.0 * m) * (0.1654 + t * (-0.0043 + t * 0.00001)) +
                sin(3.0 * m) * (0.001 + t * 0.0001) +
                cos(3.0 * m) * (0.0095 - t * 0.0003) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_2(self) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 11, 1.0)
        >>> sta2 = Saturn(epoch).station_longitude_2()
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        2018
        >>> print(m)
        9
        >>> print(round(d, 4))
        6.4175
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (68.872 + t * (-0.0007 + t * 0.00023) +
                sin(m) * (5.9399 + t * (-0.04 - t * 0.00015)) +
                cos(m) * (-0.7998 + t * (-0.0266 + t * 0.00014)) +
                sin(2.0 * m) * (0.1738 - t * 0.0032) +
                cos(2.0 * m) * (-0.0039 + t * (-0.0024 + t * 0.00001)) +
                sin(3.0 * m) * (0.0073 - t * 0.0002) +
                cos(3.0 * m) * (0.002 - t * 0.0002) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2047, 1, 1.0)
        >>> e = Saturn(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2047
        >>> print(m)
        7
        >>> print(d)
        15
        >>> print(h)
        0
        """

        # First approximation
        k = 0.03393 * (self.epoch.year() - 2003.52)
        k = round(k + 0.5) - 0.5
        jde = 2452830.12 + k * (10764.21676 - k * 0.000827)
        # Compute the epochs three months before and after
        # Compute the Sun-Saturn distance for each epoch
        sol = self._interpolate_jde(jde, delta=90.0)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(1944, 1, 1.0)
        >>> e = Saturn(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1944
        >>> print(m)
        9
        >>> print(d)
        8
        >>> print(h)
        1
        """

        # First approximation
        k = 0.03393 * (self.epoch.year() - 2003.52)
        k = round(k)
        jde = 2452830.12 + k * (10764.21676 - k * 0.000827)
        # Compute the epochs three months before and after
        # Compute the Sun-Saturn distance for each epoch
        sol = self._interpolate_jde(jde, delta=90.0)
        return Epoch(sol)

    @staticmethod
    def magnitude(sun_dist, earth_dist, delta_u, b):
        """This function computes the approximate magnitude of Saturn.

        :param sun_dist: Distance from Saturn to the Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Saturn to Earth, in Astronomical Units
        :type earth_dist: float
        :param delta_u: Difference between the Saturnicentric longitudes of the
            Sun and the Earth, measured in the plane of the ring
        :type delta_u: float, :py:class:`Angle`
        :param b: Saturnicentric latitude of the Earth refered to the plane of
            the ring, positive towards the north
        :type b: float, :py:class:`Angle`

        :returns: Saturn's magnitude
        :rtype: float

        >>> sun_dist = 9.867882
        >>> earth_dist = 10.464606
        >>> delta_u = Angle(16.442)
        >>> b = Angle(4.198)
        >>> m = Saturn.magnitude(sun_dist, earth_dist, delta_u, b)
        >>> print(m)
        1.9
        """

        # WARNING: According to Example 41.d in page 286 of Meeus book, the
        # result for the example above is 0.9 (instead of 1.9). However, after
        # carefully checking the formula implemented here, I'm sure that the
        # book has an error
        delta_u = float(delta_u)
        b = Angle(b).rad()
        m = (-8.68 + 5.0 * log10(sun_dist * earth_dist) + 0.044 * abs(delta_u)
             - 2.6 * sin(abs(b)) + 1.25 * sin(b) * sin(b))
        return round(m, 1)
