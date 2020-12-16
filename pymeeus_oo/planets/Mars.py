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


from math import sin, cos, log10

from pymeeus_oo.calculation.Angle import Angle
from pymeeus_oo.calculation.Epoch import Epoch
from pymeeus_oo.calculation.Interpolation import Interpolation
from pymeeus_oo.constellations.constellations import Constellation
from pymeeus_oo.parameters.Mars_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.planets.Earth import Earth
from pymeeus_oo.planets.Planet import Planet

"""
.. module:: Mars
   :synopsis: Class to model Mars planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Mars(Planet):
    """
    Class Mars models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    @staticmethod
    def conjunction(epoch) -> Epoch:
        """This method computes the time of the conjunction closest to the
        given epoch.

        :param epoch: Epoch close to the desired conjunction
        :type epoch: :py:class:`Epoch`

        :returns: The time when the conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(1993, 10, 1.0)
        >>> conj = Mars.conjunction(epoch)
        >>> y, m, d = conj.get_date()
        >>> print(y)
        1993
        >>> print(m)
        12
        >>> print(round(d, 4))
        27.0898
        """

        # Check that the input epoch is within valid range
        y = epoch.year()
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

    @staticmethod
    def opposition(epoch) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2729, 10, 1.0)
        >>> oppo = Mars.opposition(epoch)
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        2729
        >>> print(m)
        9
        >>> print(round(d, 4))
        9.1412
        """

        # Check that the input epoch is within valid range
        y = epoch.year()
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

    @staticmethod
    def station_longitude_1(epoch) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.
        """

        # Check that the input epoch is within valid range
        y = epoch.year()
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

    @staticmethod
    def station_longitude_2(epoch) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.
        """

        # Check that the input epoch is within valid range
        y = epoch.year()
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

    @staticmethod
    def perihelion(epoch) -> Epoch:
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

        # First approximation
        k = 0.53166 * (epoch.year() - 2001.78)
        k = round(k)  # formula for perihelion
        jde = 2452195.026 + k * (686.9957857 - k * 0.0000001187)
        # Compute the neighboring epochs half a day before and after
        sol = Mars.interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    @staticmethod
    def aphelion(epoch) -> Epoch:
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

        # First approximation
        k = 0.53166 * (epoch.year() - 2001.78)
        k = round(k + 0.5) - 0.5  # formula for aphelion
        jde = 2452195.026 + k * (686.9957857 - k * 0.0000001187)
        # Compute the neighboring epochs half a day before and after
        sol = Mars.interpolate_jde(jde, delta=0.5)
        return Epoch(sol)

    @staticmethod
    def interpolate_jde(jde: float, delta: float) -> float:
        jde_before = jde - delta
        jde_after = jde + delta
        # Compute the Sun-Mars distance for each epoch
        l, b, r_b = Mars(Epoch(jde_before)).geometric_heliocentric_position()
        l, b, r = Mars(Epoch(jde)).geometric_heliocentric_position()
        l, b, r_a = Mars(Epoch(jde_after)).geometric_heliocentric_position()
        # Call an interpolation object
        m = Interpolation([jde_before, jde, jde_after], [r_b, r, r_a])
        sol = m.minmax()
        return sol

    @staticmethod
    def magnitude(sun_dist, earth_dist, phase_angle) -> float:
        """This function computes the approximate magnitude of Mars.

        :param sun_dist: Distance from Mars to the Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Mars to Earth, in Astronomical Units
        :type earth_dist: float
        :param phase_angle: Mars phase angle
        :type phase_angle: float, :py:class:`Angle`

        :returns: Mars' magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.
        """

        # TODO: Method 'magnitude' only makes sense in the context of a 'Constellation"
        # TODO: Integrate general magnitude pattern on Constellation
        i = float(phase_angle)
        m = -1.3 + 5.0 * log10(sun_dist * earth_dist) + 0.01486 * i
        return round(m, 1)


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Mars class
    print("\n" + 35 * "*")
    print("*** Use of Mars class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    mars = Mars(epoch)
    lon, lat, r = mars.geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    mars = Mars(epoch)
    earth = Earth(epoch)
    constellation = Constellation(earth, mars)
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Mars at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    mars = Mars(epoch)
    l, a, e, i, ome, arg = mars.orbital_elements_mean_equinox()
    print_me("Mean longitude of the planet", round(l, 6))  # 288.855211
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))  # 1.52367934
    print_me("Eccentricity of the orbit", round(e, 7))  # 0.0934599
    print_me("Inclination on plane of the ecliptic", round(i, 6))  # 1.849338
    print_me("Longitude of the ascending node", round(ome, 5))  # 50.06365
    print_me("Argument of the perihelion", round(arg, 6))  # 287.202108

    print("")

    # Compute the time of the conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conj = Mars.conjunction(epoch)
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 2729/10/1
    epoch = Epoch(2729, 10, 1.0)
    oppo = Mars.opposition(epoch)
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)

    print("")

    # Compute the time of the station in longitude #1 close to 1997/3/1
    epoch = Epoch(1997, 3, 1.0)
    sta1 = Mars.station_longitude_1(epoch)
    y, m, d = sta1.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #1", date)

    # Compute the time of the station in longitude #2 close to 1997/3/1
    epoch = Epoch(1997, 3, 1.0)
    sta2 = Mars.station_longitude_2(epoch)
    y, m, d = sta2.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #2", date)

    print("")

    # Find the epoch of the Aphelion closer to 2032/1/1
    epoch = Epoch(2032, 1, 1.0)
    e = Mars.aphelion(epoch)
    y, m, d, h, mi, s = e.get_full_date()
    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'
    print_me("The Aphelion closest to 2032/1/1 will happen on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    mars = Mars(epoch)
    time, r = mars.passage_nodes()
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2019/1/15.2
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 1.4709


if __name__ == "__main__":
    main()
