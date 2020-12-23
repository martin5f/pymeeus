from pymeeus import Earth, Angle, base, Coordinates, CurveFitting, Epoch, Interpolation, Sun
from pymeeus import Jupiter
from pymeeus import Mars
from pymeeus import Mercury
from pymeeus import Minor
from pymeeus import Neptune
from pymeeus import Pluto
from pymeeus import Saturn
from pymeeus import Uranus
from pymeeus import Venus


def main():
    # Base
    Angle.main()
    base.main()
    Coordinates.main()
    CurveFitting.main()
    Epoch.main()
    Interpolation.main()

    # Sun
    Sun.main()

    # Planets
    Mercury.main()
    Venus.main()
    Earth.main()
    Mars.main()
    Jupiter.main()
    Saturn.main()
    Uranus.main()
    Neptune.main()
    Pluto.main()
    Minor.main()


if __name__ == "__main__":
    main()
