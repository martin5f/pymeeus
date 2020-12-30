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
from abc import ABC, abstractmethod

from pyplanets.core.coordinates import (
    geometric_vsop_pos, apparent_vsop_pos, orbital_elements,
    passage_nodes_elliptic, vsop_pos_terms
)
from pyplanets.core.epoch import Epoch
from pyplanets.core.interpolation import Interpolation

"""
.. module:: Planet
   :synopsis: Abstract base class for planets
   :license: GNU Lesser General Public License v3 (LGPLv3)

.. moduleauthor:: Martin Fünffinger
"""


class Planet(ABC):
    """
    Abstract base class for planets.
    """

    def __init__(self, epoch: Epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000):
        """Initializes an object of type planet for a given date (epoch) with
         its orbital parameters (ephemeredes) given by VSOP87.

        :param epoch: Epoch to which all computations of this instance refer, as an Epoch object
        :type epoch: :py:class:`Epoch`
        :raises: TypeError if epoch is of wrong type.
        """
        # First check that input value is of correct types
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")
        self.epoch: Epoch = epoch
        self.VSOP87_L = VSOP87_L
        self.VSOP87_B = VSOP87_B
        self.VSOP87_R = VSOP87_R
        self.ORBITAL_ELEM = ORBITAL_ELEM
        self.ORBITAL_ELEM_J2000 = ORBITAL_ELEM_J2000
        super().__init__()

    def geometric_heliocentric_position(self, tofk5=True):
        """This method computes the geometric heliocentric position of a  planet
        for a given epoch, using the VSOP87 theory.

        :param tofk5: Whether or not the small correction to convert to the FK5
            system will be applied or not
        :type tofk5: bool

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        """
        return geometric_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R, tofk5)

    def apparent_heliocentric_position(self, nutation=True):
        """This method computes the apparent heliocentric position of a planet
        for a given epoch, using the VSOP87 theory.

        :param nutation: Whether the nutation correction will be applied
        :type nutation: bool

        :returns: A tuple with the heliocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        """

        return apparent_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R, nutation)

    def apparent_planetcentric_sun_position(self, nutation=True):
        """Generalization of coordinate transformation applicable to all planets.
        For Planet Earth, this method returns Sun's apparent geocentric position,
        i.e. Sun.apparent_geocentric_position(...) in PyMeeus-Classic.
        Basically a simple transformation!

        :param nutation: Whether the nutation correction will be applied
        :type nutation: bool

        :returns: A tuple with Sun's geocentric longitude and latitude (as
            :py:class:`Angle` objects), and the radius vector (as a float,
            in astronomical units), in that order
        :rtype: tuple
        """
        lon, lat, r = self.apparent_heliocentric_position(nutation)
        lon = lon.to_positive() + 180.0
        lat = -lat
        return lon, lat, r

    def orbital_elements_mean_equinox(self):
        """This method computes the orbital elements of Mars for the mean
        equinox of the date for a given epoch.

        :returns: A tuple containing the following six orbital elements:
            - Mean longitude of the planet (Angle)
            - Semimajor axis of the orbit (float, astronomical units)
            - eccentricity of the orbit (float)
            - inclination on the plane of the ecliptic (Angle)
            - longitude of the ascending node (Angle)
            - argument of the perihelion (Angle)
        :rtype: tuple
        """

        return orbital_elements(self.epoch, self.ORBITAL_ELEM, self.ORBITAL_ELEM)

    def orbital_elements_j2000(self):
        """This method computes the orbital elements of Mars for the
        standard equinox J2000.0 for a given epoch.

        :returns: A tuple containing the following six orbital elements:
            - Mean longitude of the planet (Angle)
            - Semimajor axis of the orbit (float, astronomical units)
            - eccentricity of the orbit (float)
            - inclination on the plane of the ecliptic (Angle)
            - longitude of the ascending node (Angle)
            - argument of the perihelion (Angle)
        :rtype: tuple
        """

        return orbital_elements(self.epoch, self.ORBITAL_ELEM, self.ORBITAL_ELEM_J2000)

    @abstractmethod
    def aphelion(self) -> Epoch:
        """This method computes the time of Aphelion closest to the epoch of initialization.
        Subclasses are required to implement this method.

        :returns: The epoch of the desired Aphelion
        :rtype: :py:class:`Epoch`
        """
        pass

    @abstractmethod
    def perihelion(self) -> Epoch:
        """This method computes the time of Perihelion closest to the epoch of initialization.
        Subclasses are required to implement this method.

        :returns: The epoch of the desired Perihelion
        :rtype: :py:class:`Epoch`
        """
        pass

    def passage_nodes(self, ascending=True) -> (Epoch, float):
        """This function computes the time of passage by the nodes (ascending
        or descending) of the planet, nearest to the given epoch.

        :param ascending: Whether the time of passage by the ascending (True)
            or descending (False) node will be computed
        :type ascending: bool

        :returns: Tuple containing:
            - Time of passage through the node (:py:class:`Epoch`)
            - Radius vector when passing through the node (in AU, float)
        :rtype: tuple
        """

        # Get the orbital parameters
        l, a, e, i, ome, arg = self.orbital_elements_mean_equinox()
        # Compute the time of passage through perihelion
        t = self.perihelion()
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r

    def _interpolate_jde(self, jde, delta) -> float:
        """Finds the epoch corresponding to the minimum or maximum (extremum in general) distance of the planet to Sun.
        Used by perihelion() / aphelion() - methods

        :returns: The epoch of the desired Minimum / Maximum
        :rtype: float
        """

        t_b = (jde - delta - 2451545.0) / 365250.0
        t = (jde - 2451545.0) / 365250.0
        t_a = (jde + delta - 2451545.0) / 365250.0

        # Compute the Sun-Planet distance for each epoch
        r_b = vsop_pos_terms(t_b, self.VSOP87_R)
        r = vsop_pos_terms(t, self.VSOP87_R)
        r_a = vsop_pos_terms(t_a, self.VSOP87_R)

        # Call an interpolation object
        m = Interpolation([jde - delta, jde, jde + delta], [r_b, r, r_a])
        sol = m.minmax()
        return sol

