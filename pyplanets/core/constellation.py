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


import logging
import sys
from math import sin, cos, tan, atan2, sqrt, radians, acos

from pyplanets.core.angle import Angle
from pyplanets.core.coordinates import (
    nutation_longitude, true_obliquity, ecliptical2equatorial
)
from pyplanets.core.epoch import JDE2000, Epoch
from pyplanets.planets.earth import Earth
from pyplanets.planets.planet import Planet

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                    format="[%(asctime)s] %(levelname)s:%(name)s:%(message)s", datefmt="%Y-%m-%d %H:%M:%S")
_logger = logging.getLogger(__name__)

"""
.. module:: Constellation
   :synopsis: Class to model a planetary constellation as seen from Earth.
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""

class Constellation(object):

    def __init__(self, earth: Earth, planet: Planet):
        self.earth = earth  # PoV = Planet *of* View
        self.planet = planet  # PiV = Planet *in* View

    def geocentric_position(self):
        """This method computes the geocentric position of a planet (right
        ascension and declination) for the given epoch, as well as the
        elongation angle.

        :returns: A tuple containing the right ascension, the declination and
            the elongation angle as Angle objects
        :rtype: tuple
        """

        # Compute the heliocentric position of the Earth
        # Note: fixed point of time in space
        l0, b0, r0 = self.earth.geometric_heliocentric_position(tofk5=False)
        l0r, b0r = l0.rad(), b0.rad()

        x, y, z = 0.0, 0.0, 0.0
        delta_old, delta = -1.0, 5.0
        tau = 0.0
        epoch_corrected = self.earth.epoch
        count = 0

        while delta != delta_old:
            count += 1
            delta_old = delta
            tau = 0.0057755183 * delta
            # Adjust the epoch for light-time
            epoch_corrected = self.earth.epoch - tau
            self.planet = self.planet.__class__(epoch_corrected)
            # Compute the heliocentric position of the planet
            l, b, r = self.planet.geometric_heliocentric_position(tofk5=False)
            lr, br = l.rad(), b.rad()

            # Compute cartesian coordinates and distance from planet to earth
            x = r * cos(br) * cos(lr) - r0 * cos(b0r) * cos(l0r)
            y = r * cos(br) * sin(lr) - r0 * cos(b0r) * sin(l0r)
            z = r * sin(br) - r0 * sin(b0r)
            delta = sqrt(x * x + y * y + z * z)
            _logger.debug("Iteration-# / delta = {} / {}".format(count, delta))

        _logger.debug("Total iterations / tau[min] = {} / {}".format(count, tau * 24.0 * 60.0))

        # Compute longitude and latitude
        lamb = atan2(y, x)
        beta = atan2(z, sqrt(x * x + y * y))

        # Now, let's compute the aberration effect
        t = (epoch_corrected - JDE2000) / 36525
        e = 0.016708634 + t * (-0.000042037 - t * 0.0000001267)
        pie = 102.93735 + t * (1.71946 + t * 0.00046)
        pie = radians(pie)
        lon = l0 + 180.0
        lon = lon.rad()
        k = 20.49552  # The constant of aberration
        deltal1 = k * (-cos(lon - lamb) + e * cos(pie - lamb)) / cos(beta)
        deltab1 = -k * sin(beta) * (sin(lon - lamb) - e * sin(pie - lamb))
        deltal1 = Angle(0, 0, deltal1)
        deltab1 = Angle(0, 0, deltab1)

        # Correction to FK5 system
        lamb = Angle(lamb, radians=True)
        lamb = lamb.to_positive()
        beta = Angle(beta, radians=True)
        l_prime = lamb - t * (1.397 + t * 0.00031)
        deltal2 = Angle(0, 0, -0.09033)
        a = 0.03916 * (cos(l_prime.rad()) + sin(l_prime.rad()))
        a = a * tan(b.rad())
        deltal2 += Angle(0, 0, a)
        deltab2 = 0.03916 * (cos(l_prime.rad()) - sin(l_prime.rad()))
        deltab2 = Angle(0, 0, deltab2)
        # Apply the corrections
        lamb = lamb + deltal1 + deltal2
        beta = beta + deltab1 + deltab2

        # Correction for nutation
        dpsi = nutation_longitude(epoch_corrected)
        lamb += dpsi
        e = true_obliquity(epoch_corrected)
        ra, dec = ecliptical2equatorial(lamb, beta, e)

        # Let's compute the elongation angle
        """
        TODO: 
        Note #2: should correction by tau be differentiated for the Planet *and* the Sun? the time of the 
        observation (e.g. from Earth) is fixed, but the observed positions of Sun and Planet at that time are the 
        body's position some time tau before and that tau depends on the distance of Earth to the Planet and Sun 
        respectively (right?).  
        """
        earth_corrected = self.earth.__class__(epoch_corrected)
        lons, lats, rs = earth_corrected.apparent_planetcentric_sun_position()
        lambr = lamb.rad()
        lsr = lons.rad()
        betar = beta.rad()
        elon = acos(cos(betar) * cos(lambr - lsr))
        elon = Angle(elon, radians=True)
        return ra, dec, elon
