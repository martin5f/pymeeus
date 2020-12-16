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

from pymeeus_oo.calculation.Coordinates import (
    geometric_vsop_pos, apparent_vsop_pos, orbital_elements,
    passage_nodes_elliptic
)
from pymeeus_oo.calculation.Epoch import Epoch

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

    def __init__(self, epoch: Epoch, VSOP87_L, VSOP87_B, VSOP87_R, ORBITAL_ELEM, ORBITAL_ELEM_J2000):
        self.epoch: Epoch = epoch
        self.VSOP87_L = VSOP87_L
        self.VSOP87_B = VSOP87_B
        self.VSOP87_R = VSOP87_R
        self.ORBITAL_ELEM = ORBITAL_ELEM
        self.ORBITAL_ELEM_J2000 = ORBITAL_ELEM_J2000
        super().__init__()

    def set_epoch(self, epoch: Epoch):
        self.epoch = epoch

    def geometric_heliocentric_position(self, tofk5=True):
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
        return geometric_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R, tofk5)

    def apparent_heliocentric_position(self, nutation=True):
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

        return apparent_vsop_pos(self.epoch, self.VSOP87_L, self.VSOP87_B, self.VSOP87_R, nutation)

    def apparent_planetcentric_sun_position(self, nutation=True):
        """For Planet Earth, this method returns Sun's apparent geocentric position.
        Basically a simple transformation!
        """
        lon, lat, r = self.apparent_heliocentric_position(nutation)
        # lon, lat, r = Earth.apparent_heliocentric_position(epoch, nutation)
        lon = lon.to_positive() + 180.0
        lat = -lat
        return lon, lat, r

    def orbital_elements_mean_equinox(self):
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

    def orbital_elements_j2000(self):
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
        pass

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
        pass

    def passage_nodes(self, ascending=True) -> (Epoch, float):
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
        l, a, e, i, ome, arg = self.orbital_elements_mean_equinox()
        # Compute the time of passage through perihelion
        t = self.perihelion(self.epoch)
        # Get the time of passage through the node
        time, r = passage_nodes_elliptic(arg, e, a, t, ascending)
        return time, r


def main():
    pass


if __name__ == "__main__":
    main()
