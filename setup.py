from setuptools import setup, find_packages

setup(
    name='data_buttons',
    version='0.5',
    url='http://data-buttons.readthedocs.io',
    license='GPL-3.0',
    author='Thomas Williams',
    author_email='thomas.g.williams94@gmail.com',
    description='Useful Catalogue Query/Mosaicking Tools',
    packages = find_packages(),
    python_requires='>3.5',
    install_requires=['astropy','astroquery','MontagePy','photutils',
                      'drizzlepac','stwcs','crds']
)
