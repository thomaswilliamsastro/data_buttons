##########################################
Data Buttons: Quick Mosaics for You and Me
##########################################

``data_buttons`` offers a number of tools to query archives, download image data for a given source and mosaic it into a
final image.

The biggest feature of this module is the Data Button\ :sup:`TM`\, which allows you an easy way to put in an object (and
optionally an image size), and get a mosaic out. No muss, no fuss.

============
Installation
============

* Clone this repository with ``git clone https://github.com/thomaswilliamsastro/data_buttons`` into a folder of your
  choice

* Navigate to the directory that contains ``setup.py``

* Install using ``pip install -r requirements.txt -e .``

=============
Documentation
=============

.. code-block:: python

    import os
    from data_buttons import galex

    os.chdir('directory')
    galex.galex_button('M51')

This program is very simple, but documentation is available at https://data-buttons.readthedocs.io/en/latest/.

=========
Problems?
=========

Feel free to email me, thomas.g.williams94 (at) gmail.com