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
from abc import ABC, abstractmethod
from math import sin, cos, tan, acos, atan2, sqrt, radians

from pymeeus.Angle import Angle
from pymeeus.Coordinates import (
    geometric_vsop_pos, apparent_vsop_pos, orbital_elements,
    nutation_longitude, true_obliquity, ecliptical2equatorial,
    passage_nodes_elliptic
)
from pymeeus.Earth import Earth
from pymeeus.Epoch import Epoch, JDE2000
from pymeeus.Sun import Sun

"""
.. module:: Mars
   :synopsis: Class to model Mars planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Planet(ABC):
    """
    (Abstract) base class for planets.
    """

    def __init__(self, epoch, tofk5, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000):
        self.epoch = epoch
        self.tofk5 = tofk5
        self.VSOP87_L = VSOP87_L
        self.VSOP87_B = VSOP87_B
        self.VSOP87_R = VSOP87_R
        self.ORBITAL_ELEM = ORBITAL_ELEM
        self.ORBITAL_ELEM_J2000 = ORBITAL_ELEM_J2000
        super().__init__()

    def geometric_heliocentric_position(self, epoch, tofk5=True):
        """This method computes the geometric heliocentric position of planet
        Mars for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Mars position, as an Epoch object
        :type epoch: :py:class:`Epoch`
        :param tofk5: Whether or not the small correction to convert to the FK5
            system will be applied or not
        :type tofk5: bool

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.
        """
        return geometric_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R, self.tofk5)

    def apparent_heliocentric_position(self, epoch):
        """This method computes the apparent heliocentric position of planet
        Mars for a given epoch, using the VSOP87 theory.

        :param epoch: Epoch to compute Mars position, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.
        """

        return apparent_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R)

    def orbital_elements_mean_equinox(self, epoch):
        """This method computes the orbital elements of Mars for the mean
        equinox of the date for a given epoch.

        :param epoch: Epoch to compute orbital elements, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the following six orbital elements:
            - Mean longitude of the planet (Angle)
            - Semimajor axis of the orbit (float, astronomical units)
            - eccentricity of the orbit (float)
            - inclination on the plane of the ecliptic (Angle)
            - longitude of the ascending node (Angle)
            - argument of the perihelion (Angle)
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.
        """

        return orbital_elements(self.epoch, self.ORBITAL_ELEM, self.ORBITAL_ELEM)

    def orbital_elements_j2000(self, epoch):
        """This method computes the orbital elements of Mars for the
        standard equinox J2000.0 for a given epoch.

        :param epoch: Epoch to compute orbital elements, as an Epoch object
        :type epoch: :py:class:`Epoch`

        :returns: A tuple containing the following six orbital elements:
            - Mean longitude of the planet (Angle)
            - Semimajor axis of the orbit (float, astronomical units)
            - eccentricity of the orbit (float)
            - inclination on the plane of the ecliptic (Angle)
            - longitude of the ascending node (Angle)
            - argument of the perihelion (Angle)
        :rtype: tuple
        :raises: TypeError if input values are of wrong type.
        """

        return orbital_elements(self.epoch, self.ORBITAL_ELEM, self.ORBITAL_ELEM_J2000)

    def geocentric_position(self, epoch):
        """This method computes the geocentric position of Mars (right
        ascension and declination) for the given epoch, as well as the
        elongation angle.

        :param epoch_corrected: Epoch to compute geocentric position, as an Epoch object
        :type epoch_corrected: :py:class:`Epoch`

        :returns: A tuple containing the right ascension, the declination and
            the elongation angle as Angle objects
        :rtype: tuple
        :raises: TypeError if input value is of wrong type.
        """

        # Compute the heliocentric position of Mars
        l, b, r = self.geometric_heliocentric_position(self.epoch, tofk5=False)
        # Compute the heliocentric position of the Earth
        l0, b0, r0 = Earth.geometric_heliocentric_position(self.epoch, tofk5=False)
        # Convert to radians
        lr = l.rad()
        br = b.rad()
        l0r = l0.rad()
        b0r = b0.rad()
        # Compute first iteration
        x = r * cos(br) * cos(lr) - r0 * cos(b0r) * cos(l0r)
        y = r * cos(br) * sin(lr) - r0 * cos(b0r) * sin(l0r)
        z = r * sin(br) - r0 * sin(b0r)
        delta = sqrt(x * x + y * y + z * z)
        tau = 0.0057755183 * delta
        # Adjust the epoch for light-time
        # OLD: epoch -= tau
        epoch_corrected = self.epoch - tau
        # Compute again Mars coordinates with this correction
        # TODO: increase precision by more iterations (see JupiterMoons)
        l, b, r = self.geometric_heliocentric_position(epoch_corrected, tofk5=False)
        # Compute second iteration
        lr = l.rad()
        br = b.rad()
        x = r * cos(br) * cos(lr) - r0 * cos(b0r) * cos(l0r)
        y = r * cos(br) * sin(lr) - r0 * cos(b0r) * sin(l0r)
        z = r * sin(br) - r0 * sin(b0r)
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
        lons, lats, rs = Sun.apparent_geocentric_position(epoch_corrected)
        lambr = lamb.rad()
        lsr = lons.rad()
        betar = beta.rad()
        elon = acos(cos(betar) * cos(lambr - lsr))
        elon = Angle(elon, radians=True)
        return ra, dec, elon

    @abstractmethod
    def perihelion_aphelion(self, epoch, perihelion=True):
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
        """
        pass

    def passage_nodes(self, epoch, ascending=True):
        """This function computes the time of passage by the nodes (ascending
        or descending) of Mars, nearest to the given epoch.

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
        """

        # Get the orbital parameters
        l, a, e, i, ome, arg = self.orbital_elements_mean_equinox(self.epoch)
        # Compute the time of passage through perihelion
        t = self.perihelion_aphelion(self.epoch)
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r


def main():
    pass


if __name__ == "__main__":
    main()
