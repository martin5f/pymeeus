from pymeeus_oo.examples import usage_earth, usage_base, usage_angle, usage_coordinates, usage_curvefitting, \
    usage_epoch, usage_interpolation, usage_sun
from pymeeus_oo.examples import usage_jupiter
from pymeeus_oo.examples import usage_mars
from pymeeus_oo.examples import usage_mercury
from pymeeus_oo.examples import usage_minor
from pymeeus_oo.examples import usage_neptune
from pymeeus_oo.examples import usage_pluto
from pymeeus_oo.examples import usage_saturn
from pymeeus_oo.examples import usage_uranus
from pymeeus_oo.examples import usage_venus


def main():
    # Base
    usage_angle.main()
    usage_base.main()
    usage_coordinates.main()
    usage_curvefitting.main()
    usage_epoch.main()
    usage_interpolation.main()

    # Sun
    usage_sun.main()

    # Planets
    usage_mercury.main()
    usage_venus.main()
    usage_earth.main()
    usage_mars.main()
    usage_jupiter.main()
    usage_saturn.main()
    usage_uranus.main()
    usage_neptune.main()
    usage_pluto.main()
    usage_minor.main()


if __name__ == "__main__":
    main()
