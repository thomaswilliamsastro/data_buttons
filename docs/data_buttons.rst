##########################################
Data Buttons: Quick Mosaics for You and Me
##########################################

============
Introduction
============

``data_buttons`` offers a number of tools to query archives, download image data for a given source and mosaic it into a
final image.

The biggest feature of this module is the Data Button\ :sup:`TM`\, which allows you an easy way to put in an object (and
optionally an image size), and get a mosaic out. No muss, no fuss.

============
Installation
============

* Clone the ``data_buttons`` repository with ``git clone https://github.com/thomaswilliamsastro/data_buttons`` into a
  folder of your choice

* Navigate to the directory that contains ``setup.py``

* Install using ``python setup.py install``

=======================================
Example: Creating a GALEX Mosaic of M51
=======================================

Using the Data Button\ :sup:`TM`\, we'll create a GALEX mosaic of the classic poster child galaxy, M51.

.. code-block:: python

    import os
    from data_button import galex

    os.chdir('directory')
    galex.galex_button('M51')

And that's all there is to it! The program will query the necessary databases, download the data and mosaic it for you.
It'll put the final files (in this case, M51_NUV.fits and M51_FUV.fits), into the directory you specify in the
``os.chdir()`` call.