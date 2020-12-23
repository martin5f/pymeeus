from pymeeus_oo.calculation.epoch import Epoch
from pymeeus_oo.constellations.constellations import Constellation
from pymeeus_oo.planets.earth import Earth

from pymeeus_oo.planets.uranus import Uranus


def main():
    # Let's define a small helper function
    def print_me(msg, val):
        print("{}: {}".format(msg, val))

    # Let's show some uses of Uranus class
    print("\n" + 35 * "*")
    print("*** Use of Uranus class")
    print(35 * "*" + "\n")

    # Let's now compute the heliocentric position for a given epoch
    epoch = Epoch(2018, 10, 27.0)
    lon, lat, r = Uranus(epoch).geometric_heliocentric_position()
    print_me("Geometric Heliocentric Longitude", lon.to_positive())
    print_me("Geometric Heliocentric Latitude", lat)
    print_me("Radius vector", r)

    print("")

    # Compute the geocentric position for 1992/12/20:
    epoch = Epoch(1992, 12, 20.0)
    constellation = Constellation(Earth(epoch), Uranus(epoch))
    ra, dec, elon = constellation.geocentric_position()
    print_me("Right ascension", ra.ra_str(n_dec=1))
    print_me("Declination", dec.dms_str(n_dec=1))
    print_me("Elongation", elon.dms_str(n_dec=1))

    print("")

    # Print mean orbital elements for Uranus at 2065.6.24
    epoch = Epoch(2065, 6, 24.0)
    l, a, e, i, ome, arg = Uranus(epoch).orbital_elements_mean_equinox()
    print_me("Mean longitude of the planet", round(l, 6))  # 235.517526
    print_me("Semimajor axis of the orbit (UA)", round(a, 8))  # 19.21844604
    print_me("Eccentricity of the orbit", round(e, 7))  # 0.0463634
    print_me("Inclination on plane of the ecliptic", round(i, 6))  # 0.77372
    print_me("Longitude of the ascending node", round(ome, 5))  # 74.34776
    print_me("Argument of the perihelion", round(arg, 6))  # 99.630865

    print("")

    # Compute the time of the conjunction close to 1993/10/1
    epoch = Epoch(1993, 10, 1.0)
    conj = Uranus(epoch).conjunction()
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 1780/12/1
    epoch = Epoch(1780, 12, 1.0)
    oppo = Uranus(epoch).opposition()
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)

    print("")

    # Find the epoch of the Perihelion closer to 1780/1/1
    epoch = Epoch(1780, 1, 1.0)
    e = Uranus(epoch).perihelion()
    y, m, d = e.get_date()
    peri = str(y) + '/' + str(m) + '/' + str(int(d))
    print_me("The Perihelion closest to 1780/1/1 happened on", peri)

    print("")

    # Compute the time of passage through an ascending node
    epoch = Epoch(2019, 1, 1)
    time, r = Uranus(epoch).passage_nodes()
    y, m, d = time.get_date()
    d = round(d, 1)
    print("Time of passage through ascending node: {}/{}/{}".format(y, m, d))
    # 2028/8/23.2
    print("Radius vector at ascending node: {}".format(round(r, 4)))  # 19.3201


if __name__ == "__main__":
    main()
