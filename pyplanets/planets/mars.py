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
from pyplanets.parameters.mars_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Mars
   :synopsis: Class to model Mars planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Mars(Planet):
    """
    Class Mars models that planet.
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
        >>> conj = Mars(epoch).conjunction()
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1993
        >>> print(m)
        12
        >>> print(round(d, 4))
        27.0898
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mars' conjunction
        a = 2451707.414
        b = 779.936104
        m0 = 157.6047
        m1 = 48.705244
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (0.3102 + t * (-0.0001 + t * 0.00001) +
                sin(m) * (9.7273 + t * (-0.0156 + t * 0.00001)) +
                cos(m) * (-18.3195 + t * (-0.0467 + t * 0.00009)) +
                sin(2.0 * m) * (-1.6488 + t * (-0.0133 + t * 0.00001)) +
                cos(2.0 * m) * (-2.6117 + t * (-0.002 + t * 0.00004)) +
                sin(3.0 * m) * (-0.6827 + t * (-0.0026 + t * 0.00001)) +
                cos(3.0 * m) * (0.0281 + t * (0.0035 + t * 0.00001)) +
                sin(4.0 * m) * (-0.0823 + t * (0.0006 + t * 0.00001)) +
                cos(4.0 * m) * (0.1584 + t * 0.0013) +
                sin(5.0 * m) * (0.027 + t * 0.0005) +
                cos(5.0 * m) * (0.0433))
        to_return = jde0 + corr
        return Epoch(to_return)

    def opposition(self) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2729, 10, 1.0)
        >>> oppo = Mars(epoch).opposition()
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        2729
        >>> print(m)
        9
        >>> print(round(d, 4))
        9.1412
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mars' opposition
        a = 2452097.382
        b = 779.936104
        m0 = 181.9573
        m1 = 48.705244
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-0.3088 + t * t * 0.00002 +
                sin(m) * (-17.6965 + t * (0.0363 + t * 0.00005)) +
                cos(m) * (18.3131 + t * (0.0467 - t * 0.00006)) +
                sin(2.0 * m) * (-0.2162 + t * (-0.0198 - t * 0.00001)) +
                cos(2.0 * m) * (-4.5028 + t * (-0.0019 + t * 0.00007)) +
                sin(3.0 * m) * (0.8987 + t * (0.0058 - t * 0.00002)) +
                cos(3.0 * m) * (0.7666 + t * (-0.005 - t * 0.00003)) +
                sin(4.0 * m) * (-0.3636 + t * (-0.0001 + t * 0.00002)) +
                cos(4.0 * m) * (0.0402 + t * 0.0032) +
                sin(5.0 * m) * (0.0737 - t * 0.0008) +
                cos(5.0 * m) * (-0.098 - t * 0.0011))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_1(self) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mars' opposition
        a = 2452097.382
        b = 779.936104
        m0 = 181.9573
        m1 = 48.705244
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-37.079 + t * (-0.0009 + t * 0.00002) +
                sin(m) * (-20.0651 + t * (0.0228 + t * 0.00004)) +
                cos(m) * (14.5205 + t * (0.0504 - t * 0.00001)) +
                sin(2.0 * m) * (1.1737 - t * 0.0169) +
                cos(2.0 * m) * (-4.255 + t * (-0.0075 + t * 0.00008)) +
                sin(3.0 * m) * (0.4897 + t * (0.0074 - t * 0.00001)) +
                cos(3.0 * m) * (1.1151 + t * (-0.0021 - t * 0.00005)) +
                sin(4.0 * m) * (-0.3636 + t * (-0.002 + t * 0.00001)) +
                cos(4.0 * m) * (-0.1769 + t * (0.0028 + t * 0.00002)) +
                sin(5.0 * m) * (0.1437 - t * 0.0004) +
                cos(5.0 * m) * (-0.0383 - t * 0.0016))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_2(self) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Mars' opposition
        a = 2452097.382
        b = 779.936104
        m0 = 181.9573
        m1 = 48.705244
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (36.7191 + t * (0.0016 + t * 0.00003) +
                sin(m) * (-12.6163 + t * (0.0417 - t * 0.00001)) +
                cos(m) * (20.1218 + t * (0.0379 - t * 0.00006)) +
                sin(2.0 * m) * (-1.636 - t * 0.019) +
                cos(2.0 * m) * (-3.9657 + t * (0.0045 + t * 0.00007)) +
                sin(3.0 * m) * (1.1546 + t * (0.0029 - t * 0.00003)) +
                cos(3.0 * m) * (0.2888 + t * (-0.0073 - t * 0.00002)) +
                sin(4.0 * m) * (-0.3128 + t * (0.0017 + t * 0.00002)) +
                cos(4.0 * m) * (0.2513 + t * (0.0026 - t * 0.00002)) +
                sin(5.0 * m) * (-0.0021 - t * 0.0016) +
                cos(5.0 * m) * (-0.1497 - t * 0.0006))
        to_return = jde0 + corr
        return Epoch(to_return)

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2032, 1, 1.0)
        >>> e = Mars(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2032
        >>> print(m)
        10
        >>> print(d)
        24
        >>> print(h)
        22
        """

        # First approximation
        k = 0.53166 * (self.epoch.year() - 2001.78)
        k = round(k + 0.5) - 0.5  # formula for aphelion
        jde = 2452195.026 + k * (686.9957857 - k * 0.0000001187)
        # Compute the neighboring epochs half a day before and after
        sol = self._interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion (or Aphelion) closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2019, 2, 23.0)
        >>> e = Mars(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2018
        >>> print(m)
        9
        >>> print(d)
        16
        >>> print(h)
        12
        """

        # First approximation
        k = 0.53166 * (self.epoch.year() - 2001.78)
        k = round(k)  # formula for perihelion
        jde = 2452195.026 + k * (686.9957857 - k * 0.0000001187)
        # Compute the neighboring epochs half a day before and after
        sol = self._interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    def magnitude(self, sun_dist, earth_dist, phase_angle) -> float:
        """This function computes the approximate magnitude of Mars.

        :param sun_dist: Distance from Mars to the Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Mars to Earth, in Astronomical Units
        :type earth_dist: float
        :param phase_angle: Mars phase angle
        :type phase_angle: float, :py:class:`Angle`

        :returns: Mars' magnitude
        :rtype: float
        """

        # TODO: Method 'magnitude' only makes sense in the context of a 'Constellation"
        # TODO: Integrate general magnitude pattern on Constellation
        i = float(phase_angle)
        m = -1.3 + 5.0 * log10(sun_dist * earth_dist) + 0.01486 * i
        return round(m, 1)
