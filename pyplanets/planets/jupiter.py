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
from pyplanets.parameters.jupiter_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Jupiter
   :synopsis: Class to model Jupiter planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Jupiter(Planet):
    """
    Class Jupiter models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    def conjunction(self) -> Epoch:
        """This method computes the time of the conjunction closest to the
        given epoch.

        :returns: The time when the conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conj = Jupiter(epoch).conjunction()
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1993
        >>> print(m)
        10
        >>> print(round(d, 4))
        18.3341
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Jupiter's conjunction
        a = 2451671.186
        b = 398.884046
        m0 = 121.898
        m1 = 33.140229
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        aa = Angle(aa).rad()    # Convert to radians
        corr = (0.1027 + t * (0.0002 - t * 0.00009) +
                sin(m) * (-2.2637 + t * (0.0163 - t * 0.00003)) +
                cos(m) * (-6.154 + t * (-0.021 + t * 0.00008)) +
                sin(2.0 * m) * (-0.2021 + t * (-0.0017 + t * 0.00001)) +
                cos(2.0 * m) * (0.131 - t * 0.0008) +
                sin(3.0 * m) * (0.0086) +
                cos(3.0 * m) * (0.0087 + t * 0.0002) +
                sin(aa) * (0.0 + t * (0.0144 - t * 0.00008)) +
                cos(aa) * (0.3642 + t * (-0.0019 - t * 0.00029)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def opposition(self) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(-6, 9, 1.0)
        >>> oppo = Jupiter(epoch).opposition()
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        -6
        >>> print(m)
        9
        >>> print(round(d, 4))
        15.2865
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Jupiter's opposition
        a = 2451870.628
        b = 398.884046
        m0 = 318.4681
        m1 = 33.140229
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        aa = Angle(aa).rad()    # Convert to radians
        corr = (-0.1029 - t * t * 0.00009 +
                sin(m) * (-1.9658 + t * (-0.0056 + t * 0.00007)) +
                cos(m) * (6.1537 + t * (0.021 - t * 0.00006)) +
                sin(2.0 * m) * (-0.2081 - t * 0.0013) +
                cos(2.0 * m) * (-0.1116 - t * 0.001) +
                sin(3.0 * m) * (0.0074 + t * 0.0001) +
                cos(3.0 * m) * (-0.0097 - t * 0.0001) +
                sin(aa) * (0.0 + t * (0.0144 - t * 0.00008)) +
                cos(aa) * (0.3642 + t * (-0.0019 - t * 0.00029)))
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
        >>> sta1 = Jupiter(epoch).station_longitude_1()
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        2018
        >>> print(m)
        3
        >>> print(round(d, 4))
        9.1288
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Jupiter's opposition
        a = 2451870.628
        b = 398.884046
        m0 = 318.4681
        m1 = 33.140229
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        aa = Angle(aa).rad()    # Convert to radians
        corr = (-60.367 + t * (-0.0001 - t * 0.00009) +
                sin(m) * (-2.3144 + t * (-0.0124 + t * 0.00007)) +
                cos(m) * (6.7439 + t * (0.0166 - t * 0.00006)) +
                sin(2.0 * m) * (-0.2259 - t * 0.001) +
                cos(2.0 * m) * (-0.1497 - t * 0.0014) +
                sin(3.0 * m) * (0.0105 + t * 0.0001) +
                cos(3.0 * m) * (-0.0098) +
                sin(aa) * (0.0 + t * (0.0144 - t * 0.00008)) +
                cos(aa) * (0.3642 + t * (-0.0019 - t * 0.00029)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_2(self) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 11, 1.0)
        >>> sta2 = Jupiter(epoch).station_longitude_2()
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        2018
        >>> print(m)
        7
        >>> print(round(d, 4))
        10.6679
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Jupiter's opposition
        a = 2451870.628
        b = 398.884046
        m0 = 318.4681
        m1 = 33.140229
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        aa = Angle(aa).rad()    # Convert to radians
        corr = (60.3023 + t * (0.0002 - t * 0.00009) +
                sin(m) * (0.3506 + t * (-0.0034 + t * 0.00004)) +
                cos(m) * (5.3635 + t * (0.0247 - t * 0.00007)) +
                sin(2.0 * m) * (-0.1872 - t * 0.0016) +
                cos(2.0 * m) * (-0.0037 - t * 0.0005) +
                sin(3.0 * m) * (0.0012 + t * 0.0001) +
                cos(3.0 * m) * (-0.0096 - t * 0.0001) +
                sin(aa) * (0.0 + t * (0.0144 - t * 0.00008)) +
                cos(aa) * (0.3642 + t * (-0.0019 - t * 0.00029)))
        to_return = jde0 + corr
        return Epoch(to_return)

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(1981, 6, 1.0)
        >>> e = Jupiter(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1981
        >>> print(m)
        7
        >>> print(d)
        28
        >>> print(h)
        6
        """

        # First approximation
        k = 0.0843 * (self.epoch.year() - 2011.2)
        k = round(k + 0.5) - 0.5
        jde = 2455636.936 + k * (4332.897065 - k * 0.0001367)
        # Compute the epochs a month before and after
        # Compute the Sun-Jupiter distance for each epoch
        sol = self._interpolate_jde(jde, delta=30.0)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2019, 2, 23.0)
        >>> e = Jupiter(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2023
        >>> print(m)
        1
        >>> print(d)
        20
        >>> print(h)
        11
        """

        # First approximation
        k = 0.0843 * (self.epoch.year() - 2011.2)
        k = round(k)
        jde = 2455636.936 + k * (4332.897065 - k * 0.0001367)
        # Compute the epochs a month before and after
        # Compute the Sun-Jupiter distance for each epoch
        sol = self._interpolate_jde(jde, delta=30.0)
        return Epoch(sol)

    @staticmethod
    def magnitude(sun_dist, earth_dist):
        """This function computes the approximate magnitude of Jupiter.

        :param sun_dist: Distance from Jupiter to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance Jupiter to Earth, in Astronomical Units
        :type earth_dist: float

        :returns: Jupiter's magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.
        """

        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)):
            raise TypeError("Invalid input types")
        m = -8.93 + 5.0 * log10(sun_dist * earth_dist)
        return round(m, 1)
