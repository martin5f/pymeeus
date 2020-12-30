# PyPlanets

> **Library of astronomical algorithms to compute planetary ephemeris written in Python**.

PyPlanets at its current state is basically a refactored version of PyMeeus. PyMeeus itself is a Python implementation
of the astronomical algorithms described in the classical book "Astronomical Algorithms, 2nd Edition, Willmann-Bell Inc.
(1998)" by Jean Meeus. The original author of PyMeeus is Dagoberto Salazar. This fork is based on
commit `9b25a79` (https://github.com/architest/pymeeus/commit/9b25a790999251afd5d68c074e144df9f74a6dfd).

## Why refactor PyMeeus?

That's a valid question of course, especially since no (or only little) functional enhancements are made. But that is
also part of the answer. Looking for astronomical libraries to "play around", I stumbled across PyMeeus which really
gets most of the job done. On the other hand, when trying to extend PyMeeus with some features I got myself lost a bit
in the functional coding so I decided to restructure and refactor the code base to my liking (and experience)
before adding new features.

Main guiding principles for refactoring have been object-orientation and a clearer separation of concerns to ease the
handling of the modules and facilitate extensibility. As a side-effect it is now e.g. straight forward to compute e.g.
the mars-centric co-ordinates of Earth for a given point in time. Well, you might correctly object here "Why would
someone require that?", but

1. it comes at no extra cost, since it's just a side effect of object orientation making it clearer, that the relative
   ecliptical position of a planet depends on the view point and therefore is not a property of a planet and
2. maybe Elon Musk will find this useful in his new home someday (@Elon: please drop me a note then).

## PyPlanets vs. PyMeeus - what's the difference?

The most obvious difference is the file (module) structure, and an increased number of modules with a lower number of
lines of code in total (hopefully increasing maintainability). The higher number of modules is due to

- moved VSOP87 parameters to separate files (../parameters)
- moved usage examples from original sources to separate modules (../usages_pymeeus)
- added common base class for planets (Planet.py)
- separated class Ellipsoid and related methods from Earth.py to ellipsoid.py
- restructured (simplified) method `vsop_pos(...)` of module Coordinates

A lot of effort has been made by Dagoberto to create meaningful tests and documentation for the PyMeeus library in the
form of test cases, usage examples and inline documentation including doctests. These are preserved to its full extent,
but needed to be updated to the new module structure. Basically it boils down to three changes:

### Replace calls to static methods with calls to instance methods

*PyMeeus:*

```python
# Let's now compute the heliocentric position for a given epoch
epoch = Epoch(2018, 10, 27.0)
lon, lat, r = Mars.geometric_heliocentric_position(epoch)
```

*PyPlanets:*

```python
# Let's now compute the heliocentric position for a given epoch
epoch = Epoch(2018, 10, 27.0)
mars = Mars(epoch)
lon, lat, r = mars.geometric_heliocentric_position()
```

### Computation of geocentric coordinates is achieved by class Constellation

*PyMeeus:*

```python
# Compute the geocentric position for 1992/12/20:
epoch = Epoch(1992, 12, 20.0)
ra, dec, elon = Mars.geocentric_position(epoch)
```

*PyPlanets:*

```python
epoch = Epoch(1992, 12, 20.0)
mars = Mars(epoch)
earth = Earth(epoch)
constellation = Constellation(earth, mars)
ra, dec, elon = constellation.geocentric_position()
```

Especially this change seems a bit verbose, but it is planned to a) generalize the Constellation class also to be useful
for viewpoints other than earth and b) introduce a convenience pattern to provide the usual geocentric coordinates.

### Method `perihelion_aphelion(...)` split into dedicated methods

*PyMeeus:*

```python
# Find the epoch of the Aphelion closer to 2032/1/1
epoch = Epoch(2032, 1, 1.0)
e = Mars.perihelion_aphelion(epoch, perihelion=False)
```

*PyPlanets:*

```python
# Find the epoch of the Aphelion closer to 2032/1/1
epoch = Epoch(2032, 1, 1.0)
e = Mars(epoch).aphelion()
```

### Functional improvements

For precise planetary positions it is necessary to introduce light-time corrections for computed and observed positions
of planets. The calculation used in the PyMeeus library computes the correction once. While in most cases probably
sufficient, the accuracy can be enhanced by recalculating the correction factor several times until the arithmetic
precision of the machine is reached. Typically, about 4 - 5 iterations are required, not having too much impact on
processing times.

### Comparison PyMeeus / PyPlanets

To make sure refactoring didn't break obvious things, two measures were taken:

- Usage examples of PyMeeus have been outsourced to ../usages_pymeeus and ported to the slightly new API. A script
  executes the examples against the "old" and "new" API and compares the output to `stdout`.
  
- Tests from PyMeeus are preserved under ../tests_pymeeus and ported to the slightly new API.

## Roadmap

Amongst other things, it is planned to

- add more algorithms from Meeus' book e.g. including ephemeris of moons
- add some visualization options, e.g. by integrating Jupyter notebooks
- have a look at [Astropy](http://www.astropy.org/) and find ways to integrate with or benefit from

## Installation

PyPlanets will be made available for easy installation via `pip install pyplanets` in the near future. In the meantime,
just grab the code from GitHub and have fun.

## Meta

Author: Martin Fünffinger

Distributed under the GNU Lesser General Public License v3 (LGPLv3). See
``LICENSE.txt`` and ``COPYING.LESSER`` for more information.

Documentation: coming soon on readthedocs.io

GitHub: [https://github.com/martin5f/pyplanets](https://github.com/martin5f/pyplanets)

If you have Sphinx installed, you can generate your own, latest documentation going to directory 'docs' and issuing:

```sh
make html
```

Then the HTML documentation pages can be found in 'build/html'.

## Contributing

The preferred method to contribute is through forking and pull requests:

1. Fork it (<https://github.com/martin5f/pyplanets/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

Please bear in mind that PyPlanets follows the PEP8 style guide for Python code
[(PEP8)](https://www.python.org/dev/peps/pep-0008/?). We suggest you install and use a linter
like [Flake8](http://flake8.pycqa.org/en/latest/) before contributing.

Additionally, PyPlanets makes heavy use of automatic tests. As a general rule, every function or method added must have
a corresponding test in the proper place in `tests` directory.

Finally, documentation is also a big thing here. Add proper and abundant documentation to your new code. This also
includes in-line comments!!!.

## What's new

* 0.3.6
    * Add method to compute rising and setting times of the Sun
* 0.3.5
    * Add method magnitude() to planet classes
* 0.3.4
    * Add method to compute the parallax correction to Earth class
* 0.3.3
    * Add methods to compute the passage through the nodes
* 0.3.2
    * Add methods to compute the perihelion and aphelion of all planets
* 0.3.1
    * Fix errors in the elongation computation, add tests and examples of use of methods 'geocentric_position()', and
      tests and examples for 'Pluto' class
* 0.3.0
    * Added 'Pluto' class
* 0.2.11
    * Added conjunction and opposition methods for Jupiter, Saturn, Uranus and Neptune
* 0.2.10
    * Added 'geocentric_position()' method to 'Minor' class, and added conjunction and opposition methods for Mercury,
      Venus and Mars.
* 0.2.9
    * Added class 'Minor', as well as functions to compute velocity of an object and length of an orbit
* 0.2.8
    * Added methods 'geocentric_position()' to all the planets
* 0.2.7
    * Added function 'kepler_equation()' to 'Coordinates' module, and 'orbital_elements' methods to classes 'Mercury', '
      Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus' and 'Neptune'
* 0.2.6
    * Added classes 'Uranus' and 'Neptune', plus additional functions in Coordinates module'
* 0.2.5
    * Added classes 'Jupiter' and 'Saturn'
* 0.2.4
    * Minor bug fixing, added methods 'ephemeris_physical_observations()' and 'beginning_synodic_rotation()', and added
      classes 'Mercury' and 'Mars'
* 0.2.3
    * Added the complete list of VSOP87 parameters to Venus, method to compute dates of equinoxes and solstices, and the
      Equation of Time
* 0.2.2
    * Added heliocentric position method for J2000.0 (Earth) and rectangular coordinates methods (Sun)
* 0.2.1
    * Added Venus module and VSOP87-based positioning methods
* 0.2.0
    * Added Sun module
* 0.1.10
    * Added functions to compute if three objects are in a straight line, and the smallest diameter of a circle
      containing them.
* 0.1.9
    * Added functions to compute relative position angles and conjunctions
* 0.1.8
    * Added functions to compute angular separation to Coordinates module
* 0.1.7
    * Added functions to compute atmospheric refraction to Coordinates module
* 0.1.6
    * Added function 'times_rise_transit_set()' to Coordinates module
* 0.1.5
    * Added functions for parallactic angle, ecliptic points in the horizon, angle between north celestial pole and
      north pole of the ecliptic, and diurnal path vs. the horizon at time of rising or setting
* 0.1.4
    * Added several conversion functions to Coordinates module
* 0.1.3
    * Added Coordinates module
* 0.1.2
    * Added precession and proper motion methods, and changed handling of Epoch class
* 0.1.1
    * Added methods related to nutation corrections
* 0.1.0
    * Earth class added
* 0.0.9
    * Significant documentation improvements
* 0.0.8
    * Epoch class finished
* 0.0.7
    * Epoch class added
* 0.0.6
    * CurveFitting class added
* 0.0.5
    * Interpolation class added
* 0.0.4
    * Angle class finished
* 0.0.3
    * Removed unnecessary dependencies
* 0.0.2
    * Documentation improvements
* 0.0.1
    * Initial commit

