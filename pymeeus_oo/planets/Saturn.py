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
from pymeeus_oo.calculation.Coordinates import (
    passage_nodes_elliptic
)
from pymeeus_oo.calculation.Epoch import Epoch
from pymeeus_oo.calculation.Interpolation import Interpolation
from pymeeus_oo.parameters.Saturn_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.parameters.Venus_params import VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000
from pymeeus_oo.planets.Planet import Planet

"""
.. module:: Saturn
   :synopsis: Class to model Saturn planet
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Dagoberto Salazar
"""


class Saturn(Planet):
    """
    Class Saturn models that planet.
    """

    def __init__(self, epoch):
        super().__init__(epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000)

    @staticmethod
    def conjunction(epoch: Epoch) -> Epoch:
        """This method computes the time of the conjunction closest to the
        given epoch.

        :param epoch: Epoch close to the desired conjunction
        :type epoch: :py:class:`Epoch`

        :returns: The time when the conjunction happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2125, 6, 1.0)
        >>> conj = Saturn.conjunction(epoch)
        >>> y, m, d = conj.get_date()
        >>> print(y)
        2125
        >>> print(m)
        8
        >>> print(round(d, 4))
        26.4035
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's conjunction
        a = 2451681.124
        b = 378.091904
        m0 = 131.6934
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute auxiliary angles
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (0.0172 + t * (-0.0006 + t * 0.00023) +
                sin(m) * (-8.5885 + t * (0.0411 + t * 0.0002)) +
                cos(m) * (-1.147 + t * (0.0352 - t * 0.00011)) +
                sin(2.0 * m) * (0.3331 + t * (-0.0034 - t * 0.00001)) +
                cos(2.0 * m) * (0.1145 + t * (-0.0045 + t * 0.00002)) +
                sin(3.0 * m) * (-0.0169 + t * 0.0002) +
                cos(3.0 * m) * (-0.0109 + t * 0.0004) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def opposition(epoch: Epoch) -> Epoch:
        """This method computes the time of the opposition closest to the given
        epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: The time when the opposition happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(-6, 9, 1.0)
        >>> oppo = Saturn.opposition(epoch)
        >>> y, m, d = oppo.get_date()
        >>> print(y)
        -6
        >>> print(m)
        9
        >>> print(round(d, 4))
        14.3709
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (-0.0209 + t * (0.0006 + t * 0.00023) +
                sin(m) * (4.5795 + t * (-0.0312 - t * 0.00017)) +
                cos(m) * (1.1462 + t * (-0.0351 + t * 0.00011)) +
                sin(2.0 * m) * (0.0985 - t * 0.0015) +
                cos(2.0 * m) * (0.0733 + t * (-0.0031 + t * 0.00001)) +
                sin(3.0 * m) * (0.0025 - t * 0.0001) +
                cos(3.0 * m) * (0.005 - t * 0.0002) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def station_longitude_1(epoch: Epoch) -> Epoch:
        """This method computes the time of the 1st station in longitude
        (i.e. when the planet is stationary and begins to move westward -
        retrograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 1st station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 11, 1.0)
        >>> sta1 = Saturn.station_longitude_1(epoch)
        >>> y, m, d = sta1.get_date()
        >>> print(y)
        2018
        >>> print(m)
        4
        >>> print(round(d, 4))
        17.9433
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (-68.884 + t * (0.0009 + t * 0.00023) +
                sin(m) * (5.5452 + t * (-0.0279 - t * 0.0002)) +
                cos(m) * (3.0727 + t * (-0.043 + t * 0.00007)) +
                sin(2.0 * m) * (0.1101 + t * (-0.0006 - t * 0.00001)) +
                cos(2.0 * m) * (0.1654 + t * (-0.0043 + t * 0.00001)) +
                sin(3.0 * m) * (0.001 + t * 0.0001) +
                cos(3.0 * m) * (0.0095 - t * 0.0003) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def station_longitude_2(epoch: Epoch) -> Epoch:
        """This method computes the time of the 2nd station in longitude
        (i.e. when the planet is stationary and begins to move eastward -
        prograde - among the starts) closest to the given epoch.

        :param epoch: Epoch close to the desired opposition
        :type epoch: :py:class:`Epoch`

        :returns: Time when the 2nd station in longitude happens, as an Epoch
        :rtype: :py:class:`Epoch`
        :raises: TypeError if input value is of wrong type.
        :raises: ValueError if input epoch outside the -2000/4000 range.

        >>> epoch = Epoch(2018, 11, 1.0)
        >>> sta2 = Saturn.station_longitude_2(epoch)
        >>> y, m, d = sta2.get_date()
        >>> print(y)
        2018
        >>> print(m)
        9
        >>> print(round(d, 4))
        6.4175
        """

        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        # Check that the input epoch is within valid range
        y = epoch.year()
        if y < -2000.0 or y > 4000.0:
            raise ValueError("Epoch outside the -2000/4000 range")
        # Set some specific constants for Saturn's opposition
        a = 2451870.17
        b = 378.091904
        m0 = 318.0172
        m1 = 12.647487
        k = round((365.2425 * y + 1721060.0 - a) / b)
        jde0 = a + k * b
        m = m0 + k * m1
        m = Angle(m).to_positive()
        m = m.rad()
        t = (jde0 - 2451545.0) / 36525.0
        # Compute an auxiliary angle
        aa = 82.74 + 40.76 * t
        bb = 29.86 + 1181.36 * t
        cc = 14.13 + 590.68 * t
        dd = 220.02 + 1262.87 * t
        # Convert to radians
        aa = Angle(aa).rad()
        bb = Angle(bb).rad()
        cc = Angle(cc).rad()
        dd = Angle(dd).rad()
        corr = (68.872 + t * (-0.0007 + t * 0.00023) +
                sin(m) * (5.9399 + t * (-0.04 - t * 0.00015)) +
                cos(m) * (-0.7998 + t * (-0.0266 + t * 0.00014)) +
                sin(2.0 * m) * (0.1738 - t * 0.0032) +
                cos(2.0 * m) * (-0.0039 + t * (-0.0024 + t * 0.00001)) +
                sin(3.0 * m) * (0.0073 - t * 0.0002) +
                cos(3.0 * m) * (0.002 - t * 0.0002) +
                sin(aa) * (0.0 + t * (-0.0337 + t * 0.00018)) +
                cos(aa) * (-0.851 + t * (0.0044 + t * 0.00068)) +
                sin(bb) * (0.0 + t * (-0.0064 + t * 0.00004)) +
                cos(bb) * (0.2397 + t * (-0.0012 - t * 0.00008)) +
                sin(cc) * (0.0 - t * 0.001) +
                cos(cc) * (0.1245 + t * 0.0006) +
                sin(dd) * (0.0 + t * (0.0024 - t * 0.00003)) +
                cos(dd) * (0.0477 + t * (-0.0005 - t * 0.00006)))
        to_return = jde0 + corr
        return Epoch(to_return)

    @staticmethod
    def perihelion_aphelion(epoch: Epoch, perihelion=True) -> Epoch:
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

        >>> epoch = Epoch(1944, 1, 1.0)
        >>> e = Saturn.perihelion_aphelion(epoch)
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        1944
        >>> print(m)
        9
        >>> print(d)
        8
        >>> print(h)
        1
        >>> epoch = Epoch(2047, 1, 1.0)
        >>> e = Saturn.perihelion_aphelion(epoch, perihelion=False)
        >>> y, m, d, h, mi, s = e.get_full_date()
        >>> print(y)
        2047
        >>> print(m)
        7
        >>> print(d)
        15
        >>> print(h)
        0
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input value")
        # First approximation
        k = 0.03393 * (epoch.year() - 2003.52)
        if perihelion:
            k = round(k)
        else:
            k = round(k + 0.5) - 0.5
        jde = 2452830.12 + k * (10764.21676 - k * 0.000827)
        # Compute the epochs three months before and after
        jde_before = jde - 90.0
        jde_after = jde + 90.0
        # Compute the Sun-Saturn distance for each epoch
        l, b, r_b = Saturn.geometric_heliocentric_position(Epoch(jde_before))
        l, b, r = Saturn.geometric_heliocentric_position(Epoch(jde))
        l, b, r_a = Saturn.geometric_heliocentric_position(Epoch(jde_after))
        # Call an interpolation object
        m = Interpolation([jde_before, jde, jde_after], [r_b, r, r_a])
        sol = m.minmax()
        return Epoch(sol)

    @staticmethod
    def passage_nodes(epoch: Epoch, ascending=True) -> (Epoch, float):
        """This function computes the time of passage by the nodes (ascending
        or descending) of Saturn, nearest to the given epoch.

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

        >>> epoch = Epoch(2019, 1, 1)
        >>> time, r = Saturn.passage_nodes(epoch)
        >>> year, month, day = time.get_date()
        >>> print(year)
        2034
        >>> print(month)
        5
        >>> print(round(day, 1))
        30.2
        >>> print(round(r, 4))
        9.0546
        """

        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input types")
        # Get the orbital parameters
        l, a, e, i, ome, arg = Saturn.orbital_elements_mean_equinox(epoch)
        # Compute the time of passage through perihelion
        t = Saturn.perihelion_aphelion(epoch)
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r

    @staticmethod
    def magnitude(sun_dist, earth_dist, delta_u, b):
        """This function computes the approximate magnitude of Saturn.

        :param sun_dist: Distance from Saturn to the Sun, in Astronomical Units
        :type sun_dist: float
        :param earth_dist: Distance from Saturn to Earth, in Astronomical Units
        :type earth_dist: float
        :param delta_u: Difference between the Saturnicentric longitudes of the
            Sun and the Earth, measured in the plane of the ring
        :type delta_u: float, :py:class:`Angle`
        :param b: Saturnicentric latitude of the Earth refered to the plane of
            the ring, positive towards the north
        :type b: float, :py:class:`Angle`

        :returns: Saturn's magnitude
        :rtype: float
        :raises: TypeError if input values are of wrong type.

        >>> sun_dist = 9.867882
        >>> earth_dist = 10.464606
        >>> delta_u = Angle(16.442)
        >>> b = Angle(4.198)
        >>> m = Saturn.magnitude(sun_dist, earth_dist, delta_u, b)
        >>> print(m)
        1.9
        """

        # WARNING: According to Example 41.d in page 286 of Meeus book, the
        # result for the example above is 0.9 (instead of 1.9). However, after
        # carefully checking the formula implemented here, I'm sure that the
        # book has an error
        if not (isinstance(sun_dist, float) and isinstance(earth_dist, float)
                and isinstance(delta_u, (float, Angle))
                and isinstance(b, (float, Angle))):
            raise TypeError("Invalid input types")
        delta_u = float(delta_u)
        b = Angle(b).rad()
        m = (-8.68 + 5.0 * log10(sun_dist * earth_dist) + 0.044 * abs(delta_u)
             - 2.6 * sin(abs(b)) + 1.25 * sin(b) * sin(b))
        return round(m, 1)


def main():

    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Saturn class
    print("\n" + 35 * "*")
    print("*** Use of Saturn class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Saturn.geometric_heliocentric_position(epoch)
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    ra, dec, elon = Saturn.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Saturn at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Saturn.orbital_elements_mean_equinox(epoch)
    print_me("Mean longitude of the planet", round(l, 6))       # 131.196871
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))   # 9.55490779
    print_me("Eccentricity of the orbit", round(e, 7))          # 0.0553209
    print_me("Inclination on plane of the ecliptic", round(i, 6))   # 2.486426
    print_me("Longitude of the ascending node", round(ome, 5))  # 114.23974
    print_me("Argument of the perihelion", round(arg, 6))       # -19.896331

    print("")

    # Compute the time of the conjunction close to 2125/6/1
    epoch = Epoch(2125, 6, 1.0)
    conj = Saturn.conjunction(epoch)
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to -6/9/1
    epoch = Epoch(-6, 9, 1.0)
    oppo = Saturn.opposition(epoch)
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)

    print("")

    # Compute the time of the station in longitude #1 close to 2018/11/1
    epoch = Epoch(2018, 11, 1.0)
    sta1 = Saturn.station_longitude_1(epoch)
    y, m, d = sta1.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #1", date)

    # Compute the time of the station in longitude #2 close to 2018/11/1
    epoch = Epoch(2018, 11, 1.0)
    sta2 = Saturn.station_longitude_2(epoch)
    y, m, d = sta2.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #2", date)

    print("")

    # Find the epoch of the Perihelion closer to 2000/1/1
    epoch = Epoch(2000, 1, 1.0)
    e = Saturn.perihelion_aphelion(epoch)
    y, m, d, h, mi, s = e.get_full_date()
    peri = str(y) + '/' + str(m) + '/' + str(d) + ' at ' + str(h) + ' hours'
    print_me("The Perihelion closest to 2000/1/1 happened on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    time, r = Saturn.passage_nodes(epoch)
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2034/5/30.2
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 9.0546


if __name__ == "__main__":

    main()
