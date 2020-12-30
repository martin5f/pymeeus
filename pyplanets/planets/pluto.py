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


from math import sin, cos, sqrt, asin, atan2

from pyplanets.core.angle import Angle
from pyplanets.core.epoch import Epoch, JDE2000
from pyplanets.sun import Sun
from pyplanets.parameters.pluto_params import PLUTO_ARGUMENT, PLUTO_LONGITUDE, PLUTO_LATITUDE, PLUTO_RADIUS_VECTOR

"""
.. module:: Pluto
   :synopsis: Class to model Pluto minor planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Pluto(object):
    """
    Class Pluto models that minor planet.
    """

    def __init__(self, epoch: Epoch):
        self.epoch = epoch

    def geometric_heliocentric_position(self, tofk5=True):
        """This method computes the geometric heliocentric position of planet
        Pluto for a given epoch.

        Note: parameter tofk *only* introduced to allow duck-typing, not used at all here

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: ValueError if input epoch outside the 1885-2099 range.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> l, b, r = Pluto(epoch).geometric_heliocentric_position()
        >>> print(round(l, 5))
        232.74071
        >>> print(round(b, 5))
        14.58782
        >>> print(round(r, 6))
        29.711111
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < 1885.0 or y > 2099.0:
            raise ValueError("Epoch outside the 1885-2099 range")
        t = (self.epoch - JDE2000) / 36525.0
        jj = 34.35 + 3034.9057 * t
        ss = 50.08 + 1222.1138 * t
        pp = 238.96 + 144.96 * t
        # Compute the arguments
        corr_lon = 0.0
        corr_lat = 0.0
        corr_rad = 0.0
        for n, argument in enumerate(PLUTO_ARGUMENT):
            iii, jjj, kkk = argument
            alpha = Angle(iii * jj + jjj * ss + kkk * pp).to_positive()
            alpha = alpha.rad()
            sin_a = sin(alpha)
            cos_a = cos(alpha)
            a_lon, b_lon = PLUTO_LONGITUDE[n]
            corr_lon += a_lon * sin_a + b_lon * cos_a
            a_lat, b_lat = PLUTO_LATITUDE[n]
            corr_lat += a_lat * sin_a + b_lat * cos_a
            a_rad, b_rad = PLUTO_RADIUS_VECTOR[n]
            corr_rad += a_rad * sin_a + b_rad * cos_a
        # The coefficients in the tables were scaled up. Let's scale them down
        corr_lon /= 1000000.0
        corr_lat /= 1000000.0
        corr_rad /= 10000000.0
        lon = Angle(238.958116 + 144.96 * t + corr_lon)
        lat = Angle(-3.908239 + corr_lat)
        radius = 40.7241346 + corr_rad
        return lon, lat, radius

    def geocentric_position(self) -> (Angle, Angle):
        """This method computes the geocentric position of Pluto (right
        ascension and declination) for the given epoch, for the standard
        equinox *J2000.0*.
        TODO: What is the reference system for the Constellation-Compute-Scheme (all other planets)?

        :returns: A tuple containing the right ascension and the declination as
            Angle objects
        :rtype: tuple
        :raises: ValueError if input epoch outside the 1885-2099 range.

        >>> epoch = Epoch(1992, 10, 13.0)
        >>> ra, dec = Pluto(epoch).geocentric_position()
        >>> print(ra.ra_str(n_dec=1))
        15h 31' 43.7''
        >>> print(dec.dms_str(n_dec=0))
        -4d 27' 29.0''
        """

        # Check that the input epoch is within valid range
        y = self.epoch.year()
        if y < 1885.0 or y > 2099.0:
            raise ValueError("Epoch outside the 1885-2099 range")
        # Compute the heliocentric position of Pluto
        ll, b, r = Pluto(self.epoch).geometric_heliocentric_position()
        # Change angles to radians
        ll = ll.rad()
        b = b.rad()
        # Values corresponding to obliquity of ecliptic (epsilon) for J2000.0
        sine = 0.397777156
        cose = 0.917482062
        x = r * cos(ll) * cos(b)
        y = r * (sin(ll) * cos(b) * cose - sin(b) * sine)
        z = r * (sin(ll) * cos(b) * sine + sin(b) * cose)
        # Compute Sun's *J2000.0* rectacngular coordinates
        xs, ys, zs = Sun.rectangular_coordinates_j2000(self.epoch)
        # Compute auxiliary quantities
        xi = x + xs
        eta = y + ys
        zeta = z + zs
        # Compute Pluto's distance to Earth
        delta = sqrt(xi * xi + eta * eta + zeta * zeta)
        # Get the light-time difference
        tau = 0.0057755183 * delta
        # Repeat the computations using the light-time correction
        ll, b, r = Pluto(self.epoch - tau).geometric_heliocentric_position()
        # Change angles to radians
        ll = ll.rad()
        b = b.rad()
        x = r * cos(ll) * cos(b)
        y = r * (sin(ll) * cos(b) * cose - sin(b) * sine)
        z = r * (sin(ll) * cos(b) * sine + sin(b) * cose)
        # Compute auxiliary quantities
        xi = x + xs
        eta = y + ys
        zeta = z + zs
        # Compute Pluto's distance to Earth
        delta = sqrt(xi * xi + eta * eta + zeta * zeta)
        # Compute right ascension and declination
        alpha = Angle(atan2(eta, xi), radians=True)
        dec = Angle(asin(zeta / delta), radians=True)
        return alpha.to_positive(), dec
