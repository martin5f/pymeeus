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
from pyplanets.parameters.neptune_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Neptune
   :synopsis: Class to model Neptune planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Neptune(Planet):
    """
    Class Neptune models that planet.
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
        >>> conj = Neptune(epoch).conjunction()
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1994
        >>> print(m)
        1
        >>> print(round(d, 4))
        11.3057
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def opposition(self) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1846, 8, 1)
        >>> oppo = Neptune(epoch).opposition()
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        1846
        >>> print(m)
        8
        >>> print(round(d, 4))
        20.1623
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
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

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closest to the epoch of initialization.
        Subclasses are required to implement this method.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`
        """
        raise NotImplementedError

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closest to the epoch of initialization.
        Subclasses are required to implement this method.

        :returns: The epoch of the desired Perihelion
        :rtype: :py:class:`Epoch`
        """
        raise NotImplementedError

    def magnitude(self, sun_dist, earth_dist):
        """This function computes the approximate magnitude of Neptune.

        :param sun_dist: Distance from Neptune to Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance Neptune to Earth, in Astronomical Units
        :type earth_dist: float

        :returns: Neptune's magnitude
        :rtype: float
        """

        m = -7.05 + 5.0 * log10(sun_dist * earth_dist)
        return round(m, 1)
