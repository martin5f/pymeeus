from pymeeus_oo.calculation.angle import Angle
from pymeeus_oo.calculation.epoch import Epoch
from pymeeus_oo.planets.minor import Minor


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Minor class
    print("\n" + 35 * "*")
    print("*** Use of Minor class")
    print(35 * "*" + "\n")

    # Let's compute the equatorial coordinates of comet Encke
    a = 2.2091404
    e = 0.8502196
    q = a * (1.0 - e)
    i = Angle(11.94524)
    omega = Angle(334.75006)
    w = Angle(186.23352)
    t = Epoch(1990, 10, 28.54502)
    epoch = Epoch(1990, 10, 6.0)
    minor = Minor(q, e, i, omega, w, t)
    ra, dec, elong = minor.geocentric_position(epoch)
    print_me("Right ascension", ra.ra_str(n_dec=1))     # 10h 34' 13.7''
    print_me("Declination", dec.dms_str(n_dec=0))       # 19d 9' 32.0''
    print_me("Elongation", round(elong, 2))             # 40.51

    print("")

    # Now compute the heliocentric ecliptical coordinates
    lon, lat = minor.heliocentric_ecliptical_position(epoch)
    print_me("Heliocentric ecliptical longitude", lon.dms_str(n_dec=1))
    # 66d 51' 57.8''
    print_me("Heliocentric ecliptical latitude", lat.dms_str(n_dec=1))
    # 11d 56' 14.4''


if __name__ == "__main__":
    main()
