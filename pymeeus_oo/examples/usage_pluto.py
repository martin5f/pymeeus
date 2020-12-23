from pymeeus_oo.calculation.epoch import Epoch
from pymeeus_oo.constellations.constellations import Constellation
from pymeeus_oo.planets.earth import Earth

from pymeeus_oo.planets.pluto import Pluto


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Pluto class
    print("\n" + 35 * "*")
    print("*** Use of Pluto class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(1992, 10, 13.0)
    lon, lat, r = Pluto(epoch).geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/10/13:
    epoch = Epoch(1992, 10, 13.0)
    ra, dec = Pluto.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))

    constellation = Constellation(Earth(epoch), Pluto(epoch))
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

if __name__ == "__main__":
    main()
