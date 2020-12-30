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


from math import sin, cos, atan2, asin

from pyplanets.core.angle import Angle
from pyplanets.core.coordinates import (
    geometric_vsop_pos
)
from pyplanets.core.ellipsoid import Ellipsoid, WGS84
from pyplanets.core.epoch import Epoch
from pyplanets.parameters.earth_params import VSOP87_L, VSOP87_B, VSOP87_R, VSOP87_L_J2000, VSOP87_B_J2000, \
    ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pyplanets.planets.planet import Planet

"""
.. module:: Earth
   :synopsis: Class to model Earth's globe
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Earth(Planet):
    """
    Class Earth models the figure of the Earth surface and, with the help of a
    configurable reference ellipsoid, provides a set of handy method to compute
    different parameters, like the distance between two points on the surface.

    Please note that here we depart a little bit from Meeus' book because the
    Earth class uses the **World Geodetic System 1984 (WGS84)** as the default
    reference ellipsoid, instead of the International Astronomical Union 1974,
    which Meeus uses. This change is done because WGS84 is regarded as more
    modern.
    """

    def __init__(self, epoch: Epoch, ellipsoid=WGS84):
        """Earth constructor.

        It takes a reference ellipsoid as input. If not provided, the ellipsoid
        used is the WGS84 by default.

        :param ellipsoid: Reference ellipsoid to be used. WGS84 by default.

        :returns: Earth object.
        :rtype: :py:class:`Earth`
        :raises: TypeError if input value is of wrong type.
        """

        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)
        self._ellip = ellipsoid

    def __str__(self):
        """Method used when trying to print the Earth object. It essentially
        returns the corresponting '__str__()' method from the reference
        ellipsoid being used.

        :returns: Semi-major equatorial radius, flattening and angular velocity
           of the current reference ellipsoid, as a string.
        :rtype: string

        >>> e = Earth(Epoch(2020, 12, 30))
        >>> s = str(e)
        >>> v = s.split(':')
        >>> print(v[0] + '|' + str(round(float(v[1]), 14)) + '|' + v[2] )
        6378137.0|0.00335281066475|7.292115e-05
        """

        return str(self._ellip)

    def __repr__(self):
        """Method providing the 'official' string representation of the object.
        It provides a valid expression that could be used to recreate the
        object.

        :returns: As string with a valid expression to recreate the object
        :rtype: string
        """

        return "{}(ellipsoid=Ellipsoid({}, {}, {}))".format(
            self.__class__.__name__, self._ellip._a, self._ellip._f,
            self._ellip._omega
        )

    def ellipsoid(self) -> Ellipsoid:
        return self._ellip

    def geometric_heliocentric_position_j2000(self, tofk5=True):
        """This method computes the geometric heliocentric position of the
        Earth for a given epoch, using the VSOP87 theory, referred to the
        equinox J2000.0.

        :param tofk5: Whether or not the small correction to convert to the FK5
            system will be applied or not
        :type tofk5: bool

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        """

        return geometric_vsop_pos(
            self.epoch, VSOP87_L_J2000, VSOP87_B_J2000, VSOP87_R, tofk5
        )

    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closer to
        a given epoch.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(2000, 4, 1.0)
        >>> e = Earth(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2000
        >>> print(m)
        7
        >>> print(d)
        3
        >>> print(h)
        23
        >>> print(mi)
        51
        >>> epoch = Epoch(2009, 5, 1.0)
        >>> e = Earth(epoch).aphelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2009
        >>> print(m)
        7
        >>> print(d)
        4
        >>> print(h)
        1
        >>> print(mi)
        41
        """

        # First approximation
        k = 0.99997 * (self.epoch.year() - 2000.01)
        k = round(k + 0.5) - 0.5
        # Compute correction to first approximation
        jde = Earth._jde_correction_aphelion(k)
        # Compute the epochs half a day before and after
        return Epoch(self._interpolate_jde(jde, delta=0.5))

    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closer to
        a given epoch.

        :returns: The epoch of the desired Perihelion
        :rtype: :py:class:`Epoch`

        >>> epoch = Epoch(1989, 11, 20.0)
        >>> e = Earth(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1990
        >>> print(m)
        1
        >>> print(d)
        4
        >>> print(h)
        17
        >>> epoch = Epoch(2003, 3, 10.0)
        >>> e = Earth(epoch).perihelion()
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2003
        >>> print(m)
        1
        >>> print(d)
        4
        >>> print(h)
        5
        >>> print(mi)
        1
        """

        # First approximation
        k = 0.99997 * (self.epoch.year() - 2000.01)
        k = round(k)
        # Compute correction to first approximation
        jde = Earth._jde_correction_perihelion(k)
        # Compute the epochs half a day before and after
        return Epoch(self._interpolate_jde(jde, delta=0.5))

    @staticmethod
    def _jde_correction_aphelion(k: float) -> float:
        a1, a2, a3, a4, a5 = Earth._jde_correction_factors(k)
        corr = (-1.352 * sin(a1.rad()) + 0.061 * sin(a2.rad())
                + 0.062 * sin(a3.rad()) + 0.029 * sin(a4.rad())
                + 0.031 * sin(a5.rad()))
        jde = 2451547.507 + k * (365.2596358 + k * 0.0000000156)
        jde += corr
        return jde

    @staticmethod
    def _jde_correction_perihelion(k: float) -> float:
        a1, a2, a3, a4, a5 = Earth._jde_correction_factors(k)
        corr = (1.278 * sin(a1.rad()) - 0.055 * sin(a2.rad())
                - 0.091 * sin(a3.rad()) - 0.056 * sin(a4.rad())
                - 0.045 * sin(a5.rad()))
        jde = 2451547.507 + k * (365.2596358 + k * 0.0000000156)
        jde += corr
        return jde

    @staticmethod
    def _jde_correction_factors(k):
        a1 = Angle(328.41 + 132.788585 * k)
        a2 = Angle(316.13 + 584.903153 * k)
        a3 = Angle(346.20 + 450.380738 * k)
        a4 = Angle(136.95 + 659.306737 * k)
        a5 = Angle(249.52 + 329.653368 * k)
        return a1, a2, a3, a4, a5

    @staticmethod
    def parallax_correction(right_ascension: Angle, declination: Angle, latitude: Angle, distance: float,
                            hour_angle: Angle, height: float = 0.0) -> (Angle, Angle):
        """This function computes the parallaxes in right ascension and
        declination in order to obtain the topocentric values.

        :param right_ascension: Geocentric right ascension, as an
            :py:class:`Angle` object
        :type right_ascension: :py:class:`Angle`
        :param declination: Geocentric declination, as an
            :py:class:`Angle` object
        :type declination: :py:class:`Angle`
        :param latitude: Latitude of the observation point
        :type latitude: :py:class:`Angle`
        :param distance: Distance from the celestial object to the Earth, in
            Astronomical Units
        :type distance: float
        :param hour_angle: Geocentric hour angle of the celestial object, as an
            :py:class:`Angle`
        :type hour_angle: :py:class:`Angle`
        :param height: Height of observation point above sea level, in meters
        :type height: float

        :returns: Tuple containing the topocentric right ascension and
            declination
        :rtype: tuple

        >>> right_ascension = Angle(22, 38, 7.25, ra=True)
        >>> declination = Angle(-15, 46, 15.9)
        >>> latitude = Angle(33, 21, 22)
        >>> distance = 0.37276
        >>> hour_angle = Angle(288.7958)
        >>> topo_ra, topo_dec = Earth.parallax_correction(right_ascension, \
                                                          declination, \
                                                          latitude, distance, \
                                                          hour_angle)
        >>> print(topo_ra.ra_str(n_dec=2))
        22h 38' 8.54''
        >>> print(topo_dec.dms_str(n_dec=1))
        -15d 46' 30.0''
        """

        # Let's start computing the equatorial horizontal parallax
        ang = Angle(0, 0, 8.794)
        sin_pi = sin(ang.rad()) / distance
        # Also, the values related to the latitude
        rho_sinphi = WGS84.rho_sinphi(latitude, height)
        rho_cosphi = WGS84.rho_cosphi(latitude, height)
        # Now, let's compute the correction for the right ascension
        delta_a = atan2(-rho_cosphi * sin_pi * sin(hour_angle.rad()),
                        cos(declination.rad()) - rho_cosphi * sin_pi *
                        cos(hour_angle.rad()))
        delta_a = Angle(delta_a, radians=True)
        # And finally, the declination already corrected
        dec = atan2((sin(declination.rad()) - rho_sinphi * sin_pi) *
                    cos(delta_a.rad()),
                    cos(declination.rad()) - rho_cosphi * sin_pi *
                    cos(hour_angle.rad()))
        dec = Angle(dec, radians=True)
        return (right_ascension + delta_a), dec

    @staticmethod
    def parallax_ecliptical(longitude: Angle, latitude: Angle, semidiameter: Angle, obs_lat: Angle,
                            obliquity: Angle, sidereal_time: Angle, distance: float, height: float = 0.0) -> (
            Angle, Angle, Angle):
        """This function computes the topocentric coordinates of a celestial
        body (Moon or planet) directly from its geocentric values in ecliptical
        coordinates.

        :param longitude: Geocentric ecliptical longitude as an
            :py:class:`Angle`
        :type longitude: :py:class:`Angle`
        :param latitude: Geocentric ecliptical latitude as an :py:class:`Angle`
        :type latitude: :py:class:`Angle`
        :param semidiameter: Geocentric semidiameter as an :py:class:`Angle`
        :type semidiameter: :py:class:`Angle`
        :param obs_lat: Latitude of the observation point
        :type obs_lat: :py:class:`Angle`
        :param obliquity: Obliquity of the eliptic, as an :py:class:`Angle`
        :type obliquity: :py:class:`Angle`
        :param sidereal_time: Local sidereal time, as an :py:class:`Angle`
        :type sidereal_time: :py:class:`Angle`
        :param distance: Distance from the celestial object to the Earth, in
            Astronomical Units
        :type distance: float
        :param height: Height of observation point above sea level, in meters
        :type height: float

        :returns: Tuple containing the topocentric longitude, latitude and
            semidiameter
        :rtype: tuple

        >>> longitude = Angle(181, 46, 22.5)
        >>> latitude = Angle(2, 17, 26.2)
        >>> semidiameter = Angle(0, 16, 15.5)
        >>> obs_lat = Angle(50, 5, 7.8)
        >>> obliquity = Angle(23, 28, 0.8)
        >>> sidereal_time = Angle(209, 46, 7.9)
        >>> distance = 0.0024650163
        >>> topo_lon, topo_lat, topo_diam = \
                Earth.parallax_ecliptical(longitude, latitude, semidiameter, \
                                          obs_lat, obliquity, sidereal_time, \
                                          distance)
        >>> print(topo_lon.dms_str(n_dec=1))
        181d 48' 5.0''
        >>> print(topo_lat.dms_str(n_dec=1))
        1d 29' 7.1''
        >>> print(topo_diam.dms_str(n_dec=1))
        16' 25.5''
        """

        # Let's start computing the equatorial horizontal parallax
        ang = Angle(0, 0, 8.794)
        sin_pi = sin(ang.rad()) / distance
        # Also, the values related to the latitude
        rho_sinphi = WGS84.rho_sinphi(obs_lat, height)
        rho_cosphi = WGS84.rho_cosphi(obs_lat, height)
        # Let's compute some auxiliary quantities
        lonr = longitude.rad()
        latr = latitude.rad()
        semir = semidiameter.rad()
        sidr = sidereal_time.rad()
        oblr = obliquity.rad()
        n = cos(lonr) * cos(latr) - rho_cosphi * sin_pi * cos(sidr)
        # Now, compute the topocentric longitude
        topo_lon = atan2(sin(lonr) * cos(latr) -
                         sin_pi * (rho_sinphi * sin(oblr) +
                                   rho_cosphi * cos(oblr) * sin(sidr)), n)
        topo_lon = Angle(topo_lon, radians=True).to_positive()
        tlonr = topo_lon.rad()
        # Compute the topocentric latitude
        topo_lat = atan2(cos(tlonr) * (sin(latr) -
                                       sin_pi * (rho_sinphi * cos(oblr) -
                                                 rho_cosphi * sin(oblr) *
                                                 sin(sidr))), n)
        topo_lat = Angle(topo_lat, radians=True).to_positive()
        # Watch out: Latitude is only valid in the +/-90 deg range
        if abs(topo_lat) > 90.0:
            topo_lat = topo_lat - 180.0
        tlatr = topo_lat.rad()
        # And finally, let's compute the topocentric semidiameter
        topo_semi = asin((cos(tlonr) * cos(tlatr) * sin(semir)) / n)
        topo_semi = Angle(topo_semi, radians=True)
        return topo_lon, topo_lat, topo_semi
