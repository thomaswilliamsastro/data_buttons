from setuptools import setup

setup(
    name='data_buttons',
    version='0.2',
    packages=['data_buttons'],
    url='http://data-buttons.readthedocs.io',
    license='GPL-3.0',
    author='Thomas Williams',
    author_email='thomas.g.williams94@gmail.com',
    description='Useful Catalogue Query/Mosaicking Tools',
    python_requires='>3.5',
    install_requires=['astropy','numpy','astroquery','MontagePy','photutils',
                      'drizzlepac','stwcs','crds']
)
