from pymeeus_oo.calculation.epoch import Epoch
from pymeeus_oo.constellations.constellations import Constellation

from pymeeus_oo.planets.earth import Earth
from pymeeus_oo.planets.mars import Mars


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
    conj = Mars(epoch).conjunction()
    y, m, d = conj.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Conjunction date", date)

    # Compute the time of the opposition close to 2729/10/1
    epoch = Epoch(2729, 10, 1.0)
    oppo = Mars(epoch).opposition()
    y, m, d = oppo.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Opposition date", date)

    print("")

    # Compute the time of the station in longitude #1 close to 1997/3/1
    epoch = Epoch(1997, 3, 1.0)
    sta1 = Mars(epoch).station_longitude_1()
    y, m, d = sta1.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #1", date)

    # Compute the time of the station in longitude #2 close to 1997/3/1
    epoch = Epoch(1997, 3, 1.0)
    sta2 = Mars(epoch).station_longitude_2()
    y, m, d = sta2.get_date()
    d = round(d, 4)
    date = "{}/{}/{}".format(y, m, d)
    print_me("Date of station in longitude #2", date)

    print("")

    # Find the epoch of the Aphelion closer to 2032/1/1
    epoch = Epoch(2032, 1, 1.0)
    e = Mars(epoch).aphelion()
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
