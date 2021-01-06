try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

config = {
    'name': 'pyplanets',
    'version': '0.4.0',
    'description': 'PyPlanets: Object-oriented refactoring of PyMeeus, a Python library '
                   'implementing astronomical algorithms.',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'keywords': 'astronomy planets ephemeris Meeus VSOP87 module library',
    'license': 'LGPLv3',
    'author': 'Martin Fünffinger',
    'author_email': 'martinfl633@gmail.com',
    'url': 'https://github.com/martin5f/pyplanets',
    'download_url': 'https://github.com/martin5f/pyplanets',
    # 'install_requires': ['nose', 'pypandoc'],
    'packages': ['pyplanets', 'pyplanets.core', 'pyplanets.parameters', 'pyplanets.planets'],
    # 'scripts': ['example.py'],
    # 'py_modules': ['base'],
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy'
    ]
}

setup(**config)
