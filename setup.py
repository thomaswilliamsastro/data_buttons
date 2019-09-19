from setuptools import setup

setup(
    name='data_buttons',
    version='0.1',
    packages=['data_buttons'],
    url='http://data-buttons.readthedocs.io',
    license='GPL-3.0',
    author='Thomas Williams',
    author_email='thomas.g.williams94@gmail.com',
    description='Useful Catalogue Query/Mosaicking Tools',
    install_requires=['astropy','numpy','astroquery','MontagePy','photutils']
)
