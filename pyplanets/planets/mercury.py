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
from pyplanets.parameters.mercury_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Mercury
   :synopsis: Class to model Mercury planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Mercury(Planet):
    """
    Class Mercury models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    def inferior_conjunction(self) -> Epoch:
        """This method computes the time of the inferior conjunction closest to
        the given epoch.

        :returns: The time when the inferior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conjunction = Mercury(epoch).inferior_conjunction()
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1993
        >>> print(m)
        11
        >>> print(round(d, 4))
        6.1449
        >>> epoch = Epoch(1631, 10, 1.0)
        >>> conjunction = Mercury(epoch).inferior_conjunction()
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1631
        >>> print(m)
        11
        >>> print(round(d, 3))
        7.306
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def superior_conjunction(self) -> Epoch:
        """This method computes the time of the superior conjunction closest to
        the given epoch.

        :returns: The time when the superior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conjunction = Mercury(epoch).superior_conjunction()
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1993
        >>> print(m)
        8
        >>> print(round(d, 4))
        29.3301
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def western_elongation(self) -> (Epoch, Angle):
        """This method computes the time of the western elongation closest to
        the given epoch, as well as the corresponding maximum elongation angle.

        :returns: A tuple with the time when the western elongation happens, as
            an Epoch, and the maximum elongation angle, as an Angle
        :rtype: tuple
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 11, 1.0)
        >>> time, elongation = Mercury(epoch).western_elongation()
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

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def eastern_elongation(self) -> (Epoch, Angle):
        """This method computes the time of the eastern elongation closest to
        the given epoch, as well as the corresponding maximum elongation angle.

        :returns: A tuple with the time when the eastern elongation happens, as
            an Epoch, and the maximum elongation angle, as an Angle
        :rtype: tuple
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1990, 8, 1.0)
        >>> time, elongation = Mercury(epoch).eastern_elongation()
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

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def station_longitude_1(self) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :returns: Time when the 1st statin in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> sta1 = Mercury(epoch).station_longitude_1()
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        1993
        >>> print(m)
        10
        >>> print(round(d, 4))
        25.9358
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def station_longitude_2(self) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> sta2 = Mercury(epoch).station_longitude_2()
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        1993
        >>> print(m)
        11
        >>> print(round(d, 4))
        15.0724
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2000, 3, 1.0)
        >>> e = Mercury(epoch).aphelion()
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

        # First approximation
        k = 4.15201 * (self.epoch.year() - 2000.12)
        k = round(k + 0.5) - 0.5
        jde = 2451590.257 + k * 87.96934963
        # Compute the epochs half a day before and after
        sol = self._interpolate_jde(jde, 0.5)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2000, 1, 1.0)
        >>> e = Mercury(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2000
        >>> print(m)
        2
        >>> print(d)
        15
        >>> print(h)
        18
        """

        # First approximation
        k = 4.15201 * (self.epoch.year() - 2000.12)
        k = round(k)
        jde = 2451590.257 + k * 87.96934963
        # Compute the epochs half a day before and after
        sol = self._interpolate_jde(jde, 0.5)
        return Epoch(sol)

    def magnitude(self, sun_dist, earth_dist, phase_angle):
        """This function computes the approximate magnitude of Mercury.

        :param sun_dist: Distance from Mercury to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance Mercury to Earth, in Astronomical Units
        :type earth_dist: float
        :param phase_angle: Mercury phase angle
        :type phase_angle: float, :py:class:`Angle`

        :returns: Mercury's magnitude
        :rtype: float
        """

        i = float(phase_angle)
        i50 = i - 50.0
        m = (1.16 + 5.0 * log10(sun_dist * earth_dist) + 0.02838 * i50 +
             0.0001023 * i50 * i50)
        return round(m, 1)
