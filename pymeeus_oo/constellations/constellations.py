from math import sin, cos, tan, atan2, sqrt, radians, acos

from pymeeus_oo.calculation.Angle import Angle
from pymeeus_oo.calculation.Coordinates import (
    nutation_longitude, true_obliquity, ecliptical2equatorial
)
from pymeeus_oo.calculation.Epoch import JDE2000
from pymeeus_oo.planets.Earth import Earth
from pymeeus_oo.planets.Planet import Planet


class Constellation(object):

    def __init__(self, earth: Earth, planet: Planet):
        self.earth = earth  ## PoV = Planet of View
        self.planet = planet  ## PiV = Planet in View

    def geocentric_position(self):
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
        l, b, r = self.planet.geometric_heliocentric_position(tofk5=False)
        print("l, b, r = {}, {}, {}".format(l, b, r))
        # Compute the heliocentric position of the Earth
        # l0, b0, r0 = Earth.geometric_heliocentric_position(self.epoch, tofk5=False)
        # earth = Earth(self.epoch, tofk5=False)
        l0, b0, r0 = self.earth.geometric_heliocentric_position(tofk5=False)
        print("l0, b0, r0 = {}, {}, {}".format(l0, b0, r0))
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
        print("tau = {}".format(tau))
        # Adjust the epoch for light-time
        # OLD: epoch -= tau
        epoch_corrected = self.planet.epoch - tau
        print("epoch-tau = {}".format(epoch_corrected))
        # Compute again Mars coordinates with this correction
        # TODO: increase precision by more iterations (see JupiterMoons)
        self.planet.set_epoch(epoch_corrected)
        l, b, r = self.planet.geometric_heliocentric_position(tofk5=False)
        print("l, b, r = {}, {}, {}".format(l, b, r))
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
        # TODO: see below
        """
        Note #1: since Sun's apparent geocentric position is just a simple transformation from Earth's 
        heliocentric position, it should be easy, to eliminate Sun here (besides the correction effect!)
        Note #2: should correction by tau be differentiated for the Planet *and* the Sun? the time of the 
        observation (e.g. from Earth) is fixed, but the observed positions of Sun and Planet at that time are the 
        body's position some time tau before and that tau depends on the distance of Earth to the Planet and Sun 
        respectively (right?).  
        """
        self.earth.set_epoch(epoch_corrected)
        lons, lats, rs = self.earth.apparent_planetcentric_sun_position()
        # lons, lats, rs = Sun.apparent_geocentric_position(epoch_corrected)
        lambr = lamb.rad()
        lsr = lons.rad()
        betar = beta.rad()
        elon = acos(cos(betar) * cos(lambr - lsr))
        elon = Angle(elon, radians=True)
        return ra, dec, elon
