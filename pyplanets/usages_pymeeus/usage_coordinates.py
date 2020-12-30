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


from math import cos

import pyplanets.core.coordinates
from pyplanets.core.angle import Angle
from pyplanets.core.epoch import Epoch, JDE2000


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Coordinate functions
    print("\n" + 35 * "*")
    print("*** Use of Coordinate functions")
    print(35 * "*" + "\n")

    # Here follows a series of important parameters related to the angle
    # between Earth's rotation axis and the ecliptic
    e0 = pyplanets.calculation.coordinates.mean_obliquity(1987, 4, 10)
    print(
        "The mean angle between Earth rotation axis and ecliptic axis for "
        + "1987/4/10 is:"
    )
    print_me("Mean obliquity", e0.dms_str(n_dec=3))  # 23d 26' 27.407''
    epsilon = pyplanets.calculation.coordinates.true_obliquity(1987, 4, 10)
    print("'True' (instantaneous) angle between those axes for 1987/4/10 is:")
    print_me("True obliquity", epsilon.dms_str(n_dec=3))  # 23d 26' 36.849''
    epsilon = pyplanets.calculation.coordinates.true_obliquity(2018, 7, 29)
    print("'True' (instantaneous) angle between those axes for 2018/7/29 is:")
    print_me("True obliquity", epsilon.dms_str(True, 4))  # 23d 26' 7.2157''

    # The nutation effect is separated in two components: One parallel to the
    # ecliptic (nutation in longitude) and other perpendicular to the ecliptic
    # (nutation in obliquity)
    print("Nutation correction in longitude for 1987/4/10:")
    dpsi = pyplanets.calculation.coordinates.nutation_longitude(1987, 4, 10)
    print_me("Nutation in longitude", dpsi.dms_str(n_dec=3))  # 0d 0' -3.788''
    print("Nutation correction in obliquity for 1987/4/10:")
    depsilon = pyplanets.calculation.coordinates.nutation_obliquity(1987, 4, 10)  # 0d 0' 9.443''
    print_me("Nutation in obliquity", depsilon.dms_str(n_dec=3))

    print("")

    # We can compute the effects of precession on the equatorial coordinates of
    # a given star, taking also into account its proper motion

    start_epoch = JDE2000
    final_epoch = Epoch(2028, 11, 13.19)
    alpha0 = Angle(2, 44, 11.986, ra=True)
    delta0 = Angle(49, 13, 42.48)  # 2h 44' 11.986''
    print_me("Initial right ascension", alpha0.ra_str(n_dec=3))
    print_me("Initial declination", delta0.dms_str(n_dec=2))  # 49d 13' 42.48''
    pm_ra = Angle(0, 0, 0.03425, ra=True)
    pm_dec = Angle(0, 0, -0.0895)
    alpha, delta = pyplanets.calculation.coordinates.precession_equatorial(
        start_epoch, final_epoch, alpha0, delta0, pm_ra, pm_dec
    )
    print_me("Final right ascension", alpha.ra_str(n_dec=3))  # 2h 46' 11.331''
    print_me("Final declination", delta.dms_str(n_dec=2))  # 49d 20' 54.54''

    print("")

    # Something similar can also be done with the ecliptical coordinates
    start_epoch = JDE2000
    final_epoch = Epoch(-214, 6, 30.0)
    lon0 = Angle(149.48194)
    lat0 = Angle(1.76549)
    print_me("Initial ecliptical longitude", round(lon0(), 5))  # 149.48194
    print_me("Initial ecliptical latitude", round(lat0(), 5))  # 1.76549
    lon, lat = pyplanets.calculation.coordinates.precession_ecliptical(start_epoch, final_epoch, lon0, lat0)
    print_me("Final ecliptical longitude", round(lon(), 3))  # 118.704
    print_me("Final ecliptical latitude", round(lat(), 3))  # 1.615

    print("")

    # It is possible to compute with relative accuracy the proper motion of the
    # stars, taking into account their distance to Sun and relative velocity
    ra = Angle(6, 45, 8.871, ra=True)
    dec = Angle(-16.716108)
    pm_ra = Angle(0, 0, -0.03847, ra=True)
    pm_dec = Angle(0, 0, -1.2053)
    dist = 2.64
    vel = -7.6
    alpha, delta = pyplanets.calculation.coordinates.motion_in_space(ra, dec, dist, vel, pm_ra, pm_dec, -1000.0)

    print_me("Right ascension, year 2000", ra.ra_str(True, 2))
    print_me("Right ascension, year 1000", alpha.ra_str(True, 2))
    # 6h 45' 47.16''
    print_me("Declination, year 2000", dec.dms_str(True, 1))
    print_me("Declination, year 1000", delta.dms_str(True, 1))
    # -16d 22' 56.0''

    print("")

    # This module provides a series of functions to convert between equatorial,
    # ecliptical, horizontal and galactic coordinates

    ra = Angle(7, 45, 18.946, ra=True)
    dec = Angle(28, 1, 34.26)
    epsilon = Angle(23.4392911)
    lon, lat = pyplanets.calculation.coordinates.equatorial2ecliptical(ra, dec, epsilon)
    print_me("Equatorial to ecliptical. Longitude", round(lon(), 5))
    # 113.21563
    print_me("Equatorial to ecliptical. Latitude", round(lat(), 5))
    # 6.68417

    print("")

    lon = Angle(113.21563)
    lat = Angle(6.68417)
    epsilon = Angle(23.4392911)
    ra, dec = pyplanets.calculation.coordinates.ecliptical2equatorial(lon, lat, epsilon)
    print_me("Ecliptical to equatorial. Right ascension", ra.ra_str(n_dec=3))
    # 7h 45' 18.946''
    print_me("Ecliptical to equatorial. Declination", dec.dms_str(n_dec=2))
    # 28d 1' 34.26''

    print("")

    lon = Angle(77, 3, 56)
    lat = Angle(38, 55, 17)
    ra = Angle(23, 9, 16.641, ra=True)
    dec = Angle(-6, 43, 11.61)
    theta0 = Angle(8, 34, 57.0896, ra=True)
    eps = Angle(23, 26, 36.87)
    # Compute correction to convert from mean to apparent sidereal time
    delta = Angle(0, 0, ((-3.868 * cos(eps.rad())) / 15.0), ra=True)
    theta0 += delta
    h = theta0 - lon - ra
    azi, ele = pyplanets.calculation.coordinates.equatorial2horizontal(h, dec, lat)
    print_me("Equatorial to horizontal: Azimuth", round(azi, 3))  # 68.034
    print_me("Equatorial to horizontal: Elevation", round(ele, 3))  # 15.125

    print("")

    azi = Angle(68.0337)
    ele = Angle(15.1249)
    lat = Angle(38, 55, 17)
    h, dec = pyplanets.calculation.coordinates.horizontal2equatorial(azi, ele, lat)
    print_me("Horizontal to equatorial. Hour angle", round(h, 4))  # 64.3521
    print_me("Horizontal to equatorial. Declination", dec.dms_str(n_dec=0))
    # -6d 43' 12.0''

    print("")

    ra = Angle(17, 48, 59.74, ra=True)
    dec = Angle(-14, 43, 8.2)
    lon, lat = pyplanets.calculation.coordinates.equatorial2galactic(ra, dec)
    print_me("Equatorial to galactic. Longitude", round(lon, 4))  # 12.9593
    print_me("Equatorial to galactic. Latitude", round(lat, 4))  # 6.0463

    print("")

    lon = Angle(12.9593)
    lat = Angle(6.0463)
    ra, dec = pyplanets.calculation.coordinates.galactic2equatorial(lon, lat)
    print_me("Galactic to equatorial. Right ascension", ra.ra_str(n_dec=1))
    # 17h 48' 59.7''
    print_me("Galactic to equatorial. Declination", dec.dms_str(n_dec=0))
    # -14d 43' 8.0''

    print("")

    # Get the ecliptic longitudes of the two points of the ecliptic which are
    # on the horizon, as well as the angle between the ecliptic and the horizon
    sidereal_time = Angle(5.0, ra=True)
    lat = Angle(51.0)
    epsilon = Angle(23.44)
    lon1, lon2, i = pyplanets.calculation.coordinates.ecliptic_horizon(sidereal_time, lat, epsilon)
    print_me(
        "Longitude of ecliptic point #1 on the horizon", lon1.dms_str(n_dec=1)
    )  # 169d 21' 29.9''
    print_me(
        "Longitude of ecliptic point #2 on the horizon", lon2.dms_str(n_dec=1)
    )  # 349d 21' 29.9''
    print_me("Angle between the ecliptic and the horizon", round(i, 0))  # 62.0

    print("")

    # Let's compute the angle of the diurnal path of a celestial body relative
    # to the horizon at the time of rising and setting
    dec = Angle(23.44)
    lat = Angle(40.0)
    j = pyplanets.calculation.coordinates.diurnal_path_horizon(dec, lat)
    print_me(
        "Diurnal path vs. horizon angle at time of rising and setting",
        j.dms_str(n_dec=1),
    )  # 45d 31' 28.4''

    print("")

    # There is a function to compute the times (in hours of the day) of rising,
    # transit and setting of a given celestial body
    longitude = Angle(71, 5, 0.0)
    latitude = Angle(42, 20, 0.0)
    alpha1 = Angle(2, 42, 43.25, ra=True)
    delta1 = Angle(18, 2, 51.4)
    alpha2 = Angle(2, 46, 55.51, ra=True)
    delta2 = Angle(18, 26, 27.3)
    alpha3 = Angle(2, 51, 7.69, ra=True)
    delta3 = Angle(18, 49, 38.7)
    h0 = Angle(-0.5667)
    delta_t = 56.0
    theta0 = Angle(11, 50, 58.1, ra=True)
    rising, transit, setting = pyplanets.calculation.coordinates.times_rise_transit_set(
        longitude,
        latitude,
        alpha1,
        delta1,
        alpha2,
        delta2,
        alpha3,
        delta3,
        h0,
        delta_t,
        theta0,
    )

    print_me("Time of rising (hours of day)", round(rising, 4))  # 12.4238
    print_me("Time of transit (hours of day)", round(transit, 3))  # 19.675
    print_me("Time of setting (hours of day, next day)", round(setting, 3))
    # 2.911

    print("")

    # The air in the atmosphere introduces an error in the elevation due to the
    # refraction. We can compute the true (airless) elevation from the apparent
    # elevation, and viceversa
    apparent_elevation = Angle(0, 30, 0.0)
    true_elevation = pyplanets.calculation.coordinates.refraction_apparent2true(apparent_elevation)
    print_me(
        "True elevation for an apparent elevation of 30'",
        true_elevation.dms_str(n_dec=1),
    )  # 1' 14.7''

    true_elevation = Angle(0, 33, 14.76)
    apparent_elevation = pyplanets.calculation.coordinates.refraction_true2apparent(true_elevation)
    print_me(
        "Apparent elevation for a true elevation of 33' 14.76''",
        apparent_elevation.dms_str(n_dec=2),
    )  # 57' 51.96''

    print("")

    # The angular separation between two celestial objects can be easily
    # computed with the 'angular_separation()' function

    alpha1 = Angle(14, 15, 39.7, ra=True)
    delta1 = Angle(19, 10, 57.0)
    alpha2 = Angle(13, 25, 11.6, ra=True)
    delta2 = Angle(-11, 9, 41.0)
    sep_ang = pyplanets.calculation.coordinates.angular_separation(alpha1, delta1, alpha2, delta2)
    print_me(
        "Angular separation between two given celestial bodies (degrees)",
        round(sep_ang, 3),
    )  # 32.793

    print("")

    # We can compute the minimum angular separation achieved between two
    # celestial objects. For that, we must provide the positions at three
    # equidistant epochs:

    # EPOCH: Sep 13th, 1978, 0h TT:
    alpha1_1 = Angle(10, 29, 44.27, ra=True)
    delta1_1 = Angle(11, 2, 5.9)
    alpha2_1 = Angle(10, 33, 29.64, ra=True)
    delta2_1 = Angle(10, 40, 13.2)
    # EPOCH: Sep 14th, 1978, 0h TT:
    alpha1_2 = Angle(10, 36, 19.63, ra=True)
    delta1_2 = Angle(10, 29, 51.7)
    alpha2_2 = Angle(10, 33, 57.97, ra=True)
    delta2_2 = Angle(10, 37, 33.4)
    # EPOCH: Sep 15th, 1978, 0h TT:
    alpha1_3 = Angle(10, 43, 1.75, ra=True)
    delta1_3 = Angle(9, 55, 16.7)
    alpha2_3 = Angle(10, 34, 26.22, ra=True)
    delta2_3 = Angle(10, 34, 53.9)
    a = pyplanets.calculation.coordinates.minimum_angular_separation(
        alpha1_1,
        delta1_1,
        alpha1_2,
        delta1_2,
        alpha1_3,
        delta1_3,
        alpha2_1,
        delta2_1,
        alpha2_2,
        delta2_2,
        alpha2_3,
        delta2_3,
    )
    # Epoch fraction:
    print_me("Minimum angular separation, epoch fraction", round(a[0], 6))
    # -0.370726
    # NOTE: Given that 'n' is negative, and Sep 14th is the middle epoch (n=0),
    # then the minimum angular separation is achieved on Sep 13th, specifically
    # at: 1.0 - 0.370726 = 0.629274 => Sep 13.629274 = Sep 13th, 15h 6' 9''

    # Minimum angular separation:
    print_me("Minimum angular separation", a[1].dms_str(n_dec=0))  # 3' 44.0''

    print("")

    # If two objects have the same right ascension, then the relative position
    # angle between them must be 0 (or 180)
    alpha1 = Angle(14, 15, 39.7, ra=True)
    delta1 = Angle(19, 10, 57.0)
    alpha2 = Angle(14, 15, 39.7, ra=True)  # Same as alpha1
    delta2 = Angle(-11, 9, 41.0)
    pos_ang = pyplanets.calculation.coordinates.relative_position_angle(alpha1, delta1, alpha2, delta2)
    print_me("Relative position angle", round(pos_ang, 1))  # 0.0

    print("")

    # Planetary conjunctions may be computed with the appropriate function
    alpha1_1 = Angle(10, 24, 30.125, ra=True)
    delta1_1 = Angle(6, 26, 32.05)
    alpha1_2 = Angle(10, 25, 0.342, ra=True)
    delta1_2 = Angle(6, 10, 57.72)
    alpha1_3 = Angle(10, 25, 12.515, ra=True)
    delta1_3 = Angle(5, 57, 33.08)
    alpha1_4 = Angle(10, 25, 6.235, ra=True)
    delta1_4 = Angle(5, 46, 27.07)
    alpha1_5 = Angle(10, 24, 41.185, ra=True)
    delta1_5 = Angle(5, 37, 48.45)
    alpha2_1 = Angle(10, 27, 27.175, ra=True)
    delta2_1 = Angle(4, 4, 41.83)
    alpha2_2 = Angle(10, 26, 32.410, ra=True)
    delta2_2 = Angle(3, 55, 54.66)
    alpha2_3 = Angle(10, 25, 29.042, ra=True)
    delta2_3 = Angle(3, 48, 3.51)
    alpha2_4 = Angle(10, 24, 17.191, ra=True)
    delta2_4 = Angle(3, 41, 10.25)
    alpha2_5 = Angle(10, 22, 57.024, ra=True)
    delta2_5 = Angle(3, 35, 16.61)
    alpha1_list = [alpha1_1, alpha1_2, alpha1_3, alpha1_4, alpha1_5]
    delta1_list = [delta1_1, delta1_2, delta1_3, delta1_4, delta1_5]
    alpha2_list = [alpha2_1, alpha2_2, alpha2_3, alpha2_4, alpha2_5]
    delta2_list = [delta2_1, delta2_2, delta2_3, delta2_4, delta2_5]
    pc = pyplanets.calculation.coordinates.planetary_conjunction(alpha1_list, delta1_list, alpha2_list,
                                                                 delta2_list)
    print_me("Epoch fraction 'n' for planetary conjunction", round(pc[0], 5))
    # 0.23797
    print_me(
        "Difference in declination at conjunction", pc[1].dms_str(n_dec=1)
    )  # 2d 8' 21.8''

    print("")

    # A planetary conjunction with a star is a little bit simpler
    alpha_1 = Angle(15, 3, 51.937, ra=True)
    delta_1 = Angle(-8, 57, 34.51)
    alpha_2 = Angle(15, 9, 57.327, ra=True)
    delta_2 = Angle(-9, 9, 3.88)
    alpha_3 = Angle(15, 15, 37.898, ra=True)
    delta_3 = Angle(-9, 17, 37.94)
    alpha_4 = Angle(15, 20, 50.632, ra=True)
    delta_4 = Angle(-9, 23, 16.25)
    alpha_5 = Angle(15, 25, 32.695, ra=True)
    delta_5 = Angle(-9, 26, 1.01)
    alpha_star = Angle(15, 17, 0.446, ra=True)
    delta_star = Angle(-9, 22, 58.47)
    alpha_list = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5]
    delta_list = [delta_1, delta_2, delta_3, delta_4, delta_5]
    pc = pyplanets.calculation.coordinates.planet_star_conjunction(alpha_list, delta_list, alpha_star,
                                                                   delta_star)
    print_me("Epoch fraction 'n' for planetary conjunction with star",
             round(pc[0], 4))  # 0.2551
    print_me("Difference in declination with star at conjunction",
             pc[1].dms_str(n_dec=0))  # 3' 38.0''

    print("")

    # It is possible to compute when a planet and two other stars will be in a
    # straight line
    alpha_1 = Angle(7, 55, 55.36, ra=True)
    delta_1 = Angle(21, 41, 3.0)
    alpha_2 = Angle(7, 58, 22.55, ra=True)
    delta_2 = Angle(21, 35, 23.4)
    alpha_3 = Angle(8, 0, 48.99, ra=True)
    delta_3 = Angle(21, 29, 38.2)
    alpha_4 = Angle(8, 3, 14.66, ra=True)
    delta_4 = Angle(21, 23, 47.5)
    alpha_5 = Angle(8, 5, 39.54, ra=True)
    delta_5 = Angle(21, 17, 51.4)
    alpha_star1 = Angle(7, 34, 16.40, ra=True)
    delta_star1 = Angle(31, 53, 51.2)
    alpha_star2 = Angle(7, 45, 0.10, ra=True)
    delta_star2 = Angle(28, 2, 12.5)
    alpha_list = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5]
    delta_list = [delta_1, delta_2, delta_3, delta_4, delta_5]
    n = pyplanets.calculation.coordinates.planet_stars_in_line(alpha_list, delta_list, alpha_star1, delta_star1,
                                                               alpha_star2, delta_star2)
    print_me("Epoch fraction 'n' when bodies are in a straight line",
             round(n, 4))  # 0.2233

    print("")

    # The function 'straight_line()' computes if three celestial bodies are in
    # line providing the angle with which the bodies differ from a great circle
    alpha1 = Angle(5, 32, 0.40, ra=True)
    delta1 = Angle(0, -17, 56.9)
    alpha2 = Angle(5, 36, 12.81, ra=True)
    delta2 = Angle(-1, 12, 7.0)
    alpha3 = Angle(5, 40, 45.52, ra=True)
    delta3 = Angle(-1, 56, 33.3)
    psi, omega = pyplanets.calculation.coordinates.straight_line(alpha1, delta1, alpha2, delta2, alpha3, delta3)
    print_me("Angle deviation from a straight line", psi.dms_str(n_dec=0))
    # 7d 31' 1.0''
    print_me("Angular distance of central point to the straight line",
             omega.dms_str(n_dec=0))  # -5' 24.0''

    print("")

    # Let's compute the size of the smallest circle that contains three bodies
    alpha1 = Angle(12, 41, 8.63, ra=True)
    delta1 = Angle(-5, 37, 54.2)
    alpha2 = Angle(12, 52, 5.21, ra=True)
    delta2 = Angle(-4, 22, 26.2)
    alpha3 = Angle(12, 39, 28.11, ra=True)
    delta3 = Angle(-1, 50, 3.7)
    d = pyplanets.calculation.coordinates.circle_diameter(alpha1, delta1, alpha2, delta2, alpha3, delta3)
    print_me(
        "Diameter of smallest circle containing three celestial bodies",
        d.dms_str(n_dec=0),
    )  # 4d 15' 49.0''

    print("")

    # Now, let's find the apparent position of a star (Theta Persei) for a
    # given epoch
    epoch = Epoch(2028, 11, 13.19)
    alpha = Angle(2, 46, 11.331, ra=True)
    delta = Angle(49, 20, 54.54)
    sun_lon = Angle(231.328)
    app_alpha, app_delta = pyplanets.calculation.coordinates.apparent_position(epoch, alpha, delta, sun_lon)
    print_me("Apparent right ascension", app_alpha.ra_str(n_dec=2))
    # 2h 46' 14.39''
    print_me("Apparent declination", app_delta.dms_str(n_dec=2))
    # 49d 21' 7.45''

    print("")

    # Convert orbital elements of an object from one equinox to another
    epoch0 = Epoch(2358042.5305)
    epoch = Epoch(2433282.4235)
    i0 = Angle(47.122)
    arg0 = Angle(151.4486)
    lon0 = Angle(45.7481)
    i1, arg1, lon1 = pyplanets.calculation.coordinates.orbital_equinox2equinox(epoch0, epoch, i0, arg0, lon0)
    print_me("New inclination", round(i1(), 3))  # 47.138
    print_me("New argument of perihelion", round(arg1(), 4))  # 151.4782
    print_me("New longitude of ascending node", round(lon1(), 4))  # 48.6037

    print("")

    # Compute the eccentric and true anomalies using Kepler's equation
    eccentricity = 0.1
    mean_anomaly = Angle(5.0)
    e, v = pyplanets.calculation.coordinates.kepler_equation(eccentricity, mean_anomaly)
    print_me("Eccentric anomaly, Case #1", round(e(), 6))  # 5.554589
    print_me("True anomaly, Case #1", round(v(), 6))  # 6.139762
    e, v = pyplanets.calculation.coordinates.kepler_equation(0.99, Angle(0.2, radians=True))
    print_me("Eccentric anomaly, Case #2", round(e(), 8))  # 61.13444578
    print_me("True anomaly, Case #2", round(v(), 6))  # 166.311977

    print("")

    # Compute the velocity of a body in a given point of its (unperturbated
    # elliptic) orbit
    r = 1.0
    a = 17.9400782
    v = pyplanets.calculation.coordinates.velocity(r, a)
    print_me("Velocity at 1 AU", round(v, 2))  # 41.53

    # Compute the velocity at perihelion
    e = 0.96727426
    vp = pyplanets.calculation.coordinates.velocity_perihelion(e, a)
    print_me("Velocity at perihelion", round(vp, 2))  # 54.52

    # Compute the velocity at aphelion
    va = pyplanets.calculation.coordinates.velocity_aphelion(e, a)
    print_me("Velocity at aphelion", round(va, 2))  # 0.91

    # Calculate the length of the orbit
    length = pyplanets.calculation.coordinates.length_orbit(e, a)
    print_me("Length of the orbit (AU)", round(length, 2))  # 77.06

    print("")

    # Passage through the nodes of an elliptic orbit
    omega = Angle(111.84644)
    e = 0.96727426
    a = 17.9400782
    t = Epoch(1986, 2, 9.45891)
    time, r = pyplanets.calculation.coordinates.passage_nodes_elliptic(omega, e, a, t)
    y, m, d = time.get_date()
    d = round(d, 2)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 1985/11/9.16
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 1.8045

    # Passage through the nodes of a parabolic orbit
    omega = Angle(154.9103)
    q = 1.324502
    t = Epoch(1989, 8, 20.291)
    time, r = pyplanets.calculation.coordinates.passage_nodes_parabolic(omega, q, t, ascending=False)
    y, m, d = time.get_date()
    d = round(d, 2)
    print("Time of passage through descending node: {}/{}/{}".format(y, m, d))
    # 1989/9/17.64
    print("Radius vector at descending node: {}".format(round(r, 4)))  # 1.3901

    print("")

    # Compute the phase angle
    sun_dist = 0.724604
    earth_dist = 0.910947
    sun_earth_dist = 0.983824
    angle = pyplanets.calculation.coordinates.phase_angle(sun_dist, earth_dist, sun_earth_dist)
    print_me("Phase angle", round(angle, 2))  # 72.96
    # Now, let's compute the illuminated fraction of the disk
    k = pyplanets.calculation.coordinates.illuminated_fraction(sun_dist, earth_dist, sun_earth_dist)
    print_me("Illuminated fraction of planet disk", round(k, 3))  # 0.647


if __name__ == "__main__":
    main()
