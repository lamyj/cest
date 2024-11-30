CEST Toolbox
============

Installation
------------

This toolbox is based on Python packages (`doit`_, `nibabel`_, `numpy`_, `scipy`_, `spire`_) and additional software (`Dicomifier`_).

The recommended way is to install the CEST toolbox through an Anaconda-compatible environment, e.g. `Miniforge <https://github.com/conda-forge/miniforge>`_. With this, the dependencies can be installed with  ``conda install -c conda-forge dicomifier nibabel numpy scipy spire-pipeline``.

Usage
-----

The toolbox includes a :doc:`simulation framework<simulation>` as well as :doc:`processing <processing>` functions (e.g. WASSR or MTR) and visualization advice.

.. image:: ./media/z_spectra.png

.. toctree::
   :caption: Contents:
   :hidden:
   
   simulation.rst
   processing.rst
   api/index.rst

.. _Dicomifier: https://dicomifier.readthedocs.io
.. _doit: https://pydoit.org
.. _nibabel: https://nipy.org/nibabel
.. _numpy: https://numpy.org
.. _scipy: https://scipy.org
.. _spire: https://github.com/lamyj/spire
