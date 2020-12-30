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
from pyplanets.parameters.uranus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Uranus
   :synopsis: Class to model Uranus planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Uranus(Planet):
    """
    Class Uranus models that planet.
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
        >>> conj = Uranus(epoch).conjunction()
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1994
        >>> print(m)
        1
        >>> print(round(d, 4))
        12.7365
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def opposition(self) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1780, 12, 1.0)
        >>> oppo = Uranus(epoch).opposition()
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        1780
        >>> print(m)
        12
        >>> print(round(d, 4))
        17.5998
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        .. note:: The solution provided by this method may have several days of
            error.

        >>> epoch = Epoch(2090, 1, 1.0)
        >>> e = Uranus(epoch).aphelion()
        >>> y, m, d = e.get_date()
        >>> print(y)
        2092
        >>> print(m)
        11
        >>> print(int(d))
        22
        """

        # First approximation
        k = 0.0119 * (self.epoch.year() - 2051.1)
        k = round(k + 0.5) - 0.5
        jde = 2470213.5 + k * (30694.8767 - k * 0.00541)
        # Compute the epochs 1 year before and after
        # Compute the Sun-Uranus distance for each epoch
        sol = self._interpolate_jde(jde, delta=360.0)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion (or Aphelion)
        :rtype: :py:class:`Epoch`

        .. note:: The solution provided by this method may have several days of
            error.

        >>> epoch = Epoch(1880, 1, 1.0)
        >>> e = Uranus(epoch).perihelion()
        >>> y, m, d = e.get_date()
        >>> print(y)
        1882
        >>> print(m)
        3
        >>> print(int(d))
        18
        """

        # First approximation
        k = 0.0119 * (self.epoch.year() - 2051.1)
        k = round(k)
        jde = 2470213.5 + k * (30694.8767 - k * 0.00541)
        # Compute the epochs 1 year before and after
        # Compute the Sun-Uranus distance for each epoch
        sol = self._interpolate_jde(jde, delta=360.0)
        return Epoch(sol)

    @staticmethod
    def magnitude(sun_dist, earth_dist):
        """This function computes the approximate magnitude of Uranus.

        :param sun_dist: Distance from Uranus to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Uranus to Earth, in Astronomical Units
        :type earth_dist: float

        :returns: Uranus's magnitude
        :rtype: float
        """

        m = -6.85 + 5.0 * log10(sun_dist * earth_dist)
        return round(m, 1)
