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


from math import sin, cos, sqrt, log10

from pyplanets.core.angle import Angle
from pyplanets.core.epoch import Epoch
from pyplanets.parameters.venus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Venus
   :synopsis: Class to model Venus planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Venus(Planet):
    """
    Class Venus models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    def inferior_conjunction(self) -> Epoch:
        """This method computes the time of the inferior conjunction closest to
        the given epoch.

        :returns: The time when the inferior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1882, 12, 1.0)
        >>> conjunction = Venus(epoch).inferior_conjunction()
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1882
        >>> print(m)
        12
        >>> print(round(d, 1))
        6.7
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' inferior conjunction
        a = 2451996.706
        b = 583.921361
        m0 = 82.7311
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-0.0096 + t * (0.0002 - t * 0.00001) +
                sin(m) * (2.0009 + t * (-0.0033 - t * 0.00001)) +
                cos(m) * (0.598 + t * (-0.0104 + t * 0.00001)) +
                sin(2.0 * m) * (0.0967 + t * (-0.0018 - t * 0.00003)) +
                cos(2.0 * m) * (0.0913 + t * (0.0009 - t * 0.00002)) +
                sin(3.0 * m) * (0.0046 - t * 0.0002) +
                cos(3.0 * m) * (0.0079 + t * 0.0001))
        to_return = jde0 + corr
        return Epoch(to_return)

    def superior_conjunction(self) -> Epoch:
        """This method computes the time of the superior conjunction closest to
        the given epoch.

        :returns: The time when the superior conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conjunction = Venus(epoch).superior_conjunction()
        >>> y, m, d = conjunction.get_date()
        >>> print(y)
        1994
        >>> print(m)
        1
        >>> print(round(d, 2))
        17.05
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' superior conjunction
        a = 2451704.746
        b = 583.921361
        m0 = 154.9745
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (0.0099 + t * (-0.0002 - t * 0.00001) +
                sin(m) * (4.1991 + t * (-0.0121 - t * 0.00003)) +
                cos(m) * (-0.6095 + t * (0.0102 - t * 0.00002)) +
                sin(2.0 * m) * (0.25 + t * (-0.0028 - t * 0.00003)) +
                cos(2.0 * m) * (0.0063 + t * (0.0025 - t * 0.00002)) +
                sin(3.0 * m) * (0.0232 + t * (-0.0005 - t * 0.00001)) +
                cos(3.0 * m) * (0.0031 + t * 0.0004))
        to_return = jde0 + corr
        return Epoch(to_return)

    def western_elongation(self) -> (Epoch, Angle):
        """This method computes the time of the western elongation closest to
        the given epoch, as well as the corresponding maximum elongation angle.

        :returns: A tuple with the time when the western elongation happens, as
            an Epoch, and the maximum elongation angle, as an Angle
        :rtype: tuple
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2019, 1, 1.0)
        >>> time, elongation = Venus(epoch).western_elongation()
        >>> y, m, d = time.get_date()
        >>> print(y)
        2019
        >>> print(m)
        1
        >>> print(round(d, 4))
        6.1895
        >>> print(round(elongation, 4))
        46.9571
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' inferior conjunction
        a = 2451996.706
        b = 583.921361
        m0 = 82.7311
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (70.7462 - t * t * 0.00001 +
                sin(m) * (1.1218 + t * (-0.0025 - t * 0.00001)) +
                cos(m) * (0.4538 - t * 0.0066) +
                sin(2.0 * m) * (0.132 + t * (0.002 - t * 0.00003)) +
                cos(2.0 * m) * (-0.0702 + t * (0.0022 + t * 0.00004)) +
                sin(3.0 * m) * (0.0062 - t * 0.0001) +
                cos(3.0 * m) * (0.0015 - t * t * 0.00001))
        elon = (46.3245 +
                sin(m) * (-0.5366 + t * (-0.0003 + t * 0.00001)) +
                cos(m) * (0.3097 + t * (0.0016 - t * 0.00001)) +
                sin(2.0 * m) * (-0.0163) +
                cos(2.0 * m) * (-0.0075 + t * 0.0001))
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

        >>> epoch = Epoch(2019, 10, 1.0)
        >>> time, elongation = Venus(epoch).eastern_elongation()
        >>> y, m, d = time.get_date()
        >>> print(y)
        2020
        >>> print(m)
        3
        >>> print(round(d, 4))
        24.9179
        >>> print(round(elongation, 4))
        46.078
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' inferior conjunction
        a = 2451996.706
        b = 583.921361
        m0 = 82.7311
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-70.76 + t * (0.0002 - t * 0.00001) +
                sin(m) * (1.0282 + t * (-0.001 - t * 0.00001)) +
                cos(m) * (0.2761 - t * 0.006) +
                sin(2.0 * m) * (-0.0438 + t * (-0.0023 + t * 0.00002)) +
                cos(2.0 * m) * (0.166 + t * (-0.0037 - t * 0.00004)) +
                sin(3.0 * m) * (0.0036 + t * 0.0001) +
                cos(3.0 * m) * (-0.0011 + t * t * 0.00001))
        elon = (46.3173 + t * 0.0001 +
                sin(m) * (0.6916 - t * 0.0024) +
                cos(m) * (0.6676 - t * 0.0045) +
                sin(2.0 * m) * (0.0309 - t * 0.0002) +
                cos(2.0 * m) * (0.0036 - t * 0.0001))
        elon = Angle(elon).to_positive()
        to_return = jde0 + corr
        return Epoch(to_return), elon

    def station_longitude_1(self) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 12, 1.0)
        >>> sta1 = Venus(epoch).station_longitude_1()
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        2018
        >>> print(m)
        10
        >>> print(round(d, 4))
        5.7908
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' inferior conjunction
        a = 2451996.706
        b = 583.921361
        m0 = 82.7311
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (-21.0672 + t * (0.0002 - t * 0.00001) +
                sin(m) * (1.9396 + t * (-0.0029 - t * 0.00001)) +
                cos(m) * (1.0727 - t * 0.0102) +
                sin(2.0 * m) * (0.0404 + t * (-0.0023 - t * 0.00001)) +
                cos(2.0 * m) * (0.1305 + t * (-0.0004 - t * 0.00003)) +
                sin(3.0 * m) * (-0.0007 - t * 0.0002) +
                cos(3.0 * m) * (0.0098))
        to_return = jde0 + corr
        return Epoch(to_return)

    def station_longitude_2(self) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 12, 1.0)
        >>> sta2 = Venus(epoch).station_longitude_2()
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        2018
        >>> print(m)
        11
        >>> print(round(d, 4))
        16.439
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Venus' inferior conjunction
        a = 2451996.706
        b = 583.921361
        m0 = 82.7311
        m1 = 215.513058
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        corr = (21.0623 - t * t * 0.00001 +
                sin(m) * (1.9913 + t * (-0.004 - t * 0.00001)) +
                cos(m) * (-0.0407 - t * 0.0077) +
                sin(2.0 * m) * (0.1351 + t * (-0.0009 - t * 0.00004)) +
                cos(2.0 * m) * (0.0303 + t * 0.0019) +
                sin(3.0 * m) * (0.0089 - t * 0.0002) +
                cos(3.0 * m) * (0.0043 + t * 0.0001))
        to_return = jde0 + corr
        return Epoch(to_return)

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(1979, 2, 1.0)
        >>> e = Venus(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1979
        >>> print(m)
        4
        >>> print(d)
        22
        >>> print(h)
        12
        """

        # First approximation
        k = 1.62549 * (self.epoch.year() - 2000.53)
        k = round(k + 0.5) - 0.5
        jde = 2451738.233 + k * (224.7008188 - k * 0.0000000327)
        # Compute the neighboring epochs half a day before and after
        sol = self._interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(1978, 10, 15.0)
        >>> e = Venus(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1978
        >>> print(m)
        12
        >>> print(d)
        31
        >>> print(h)
        4
        """

        # First approximation
        k = 1.62549 * (self.epoch.year() - 2000.53)
        k = round(k)
        jde = 2451738.233 + k * (224.7008188 - k * 0.0000000327)
        # Compute the neighboring epochs half a day before and after
        sol = self._interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    @staticmethod
    def illuminated_fraction(epoch: Epoch):
        """This function computes an approximation of the illuminated fraction
        of Venus disk, as seen from Earth.

        :param epoch: Epoch to compute the illuminated fraction
        :type epoch: :py:class:`Epoch`

        :returns: Illuminated fraction of Venus disk
        :rtype: float
        :raises: TypeError if input values are of wrong type.

        >>> epoch = Epoch(1992, 12, 20)
        >>> k = Venus.illuminated_fraction(epoch)
        >>> print(round(k, 2))
        0.64
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        t = (epoch.jde() - 2451545.0) / 36525.0
        v = Angle(261.51 + 22518.443 * t)
        m = Angle(177.53 + 35999.050 * t)
        n = Angle(50.42 + 58517.811 * t)
        w = Angle(v + 1.91 * sin(m.rad()) + 0.78 * sin(n.rad()))
        delta2 = abs(1.52321 + 1.44666 * cos(w.rad()))
        delta = sqrt(delta2)
        k = ((0.72333 + delta) * (0.72333 + delta) - 1.0) / (2.89332 * delta)
        return k

    @staticmethod
    def magnitude(sun_dist, earth_dist, phase_angle):
        """This function computes the approximate magnitude of Venus.

        :param sun_dist: Distance from Venus to the Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Venus to Earth, in Astronomical Units
        :type earth_dist: float
        :param phase_angle: Venus phase angle
        :type phase_angle: float, :py:class:`Angle`

        :returns: Venus' magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.

        >>> sun_dist = 0.724604
        >>> earth_dist = 0.910947
        >>> phase_angle = Angle(72.96)
        >>> m = Venus.magnitude(sun_dist, earth_dist, phase_angle)
        >>> print(m)
        -3.8
        """

        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)
                and isinstance(phase_angle, (float, Angle))):
            raise TypeError("Invalid input types")
        i = float(phase_angle)
        m = (-4.0 + 5.0 * log10(sun_dist * earth_dist) + 0.01322 * i
             + 0.0000004247 * i * i * i)
        return round(m, 1)
