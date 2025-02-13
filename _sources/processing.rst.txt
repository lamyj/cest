Processing
==========

There are two main ways to use the post-processing features this toolbox: through function calls, or through tasks in a `spire`_ pipeline. The former approach works well for experimenting and when processing a limited quantity of data while the latter approach is particularly suited for reproducible analysis on a larger database.

Functions
---------

Functions work on files: they read and write data from the disk, and thus have no return value. Images, including multi-volume images of Z-spectra, are expected to be stored in `well-known format <https://nipy.org/nibabel/api.html>`_. Meta-data are expected in the JSON format used by `Dicomifier`_.

For example, a shift map using the `WASSR`_ method can be computed using the following code sample. It will load data from a Z-spectrum map and its associated meta-data, and store the result in ``deta_ppm.nii.gz``.

.. code:: python
   
   import cest
   
   cest.wassr("z_spectrum.nii.gz", "z_spectrum.json", "delta_ppm.nii.gz")

The documentation of individual functions is described in the :doc:`API documentation <api/functions>`.

Tasks
-----

By linking tasks in a `spire`_ pipeline, the depencies between processing steps are handled automatically, improving the reproducibility of your pipeline. The example below assumes that individual exams are stored as subdirectories of ``/storage/my_study``, and, for each exam, creates a shift map in PPM, shifts the Z-spectra and resamples it in the [-5, +5] PPM range.

.. code:: python
   
   import pathlib
   
   import cest.tasks
   import numpy
   
   root = pathlib.Path("/storage/my_study")
   for exam in root.iterdir():
       shift_map = cest.tasks.WASSR(
           exam/"water.nii.gz", exam/"water.json", exam/"delta_ppm.nii.gz")
       shifted = cest.tasks.ShiftSpectrum(
           exam/"glutamate.nii.gz", exam/"glutamate.json", shift_map.targets[0],
           exam/"glutamate_shifted.nii.gz")
       refined = cest.tasks.Refine(
           shifted.targets[0], shifted.file_dep[1], numpy.linspace(-5, 5, 501),
           derived/"Glutamate_refined.nii.gz")

The documentation of tasks is available in the :doc:`API documentation <api/tasks>`.

Visualization
-------------

`FSLeyes`_ can load NIfTI images and display the Z-spectrum of a single voxel. The provided `layout file <_static/cest.txt>`_ can be `added to the pre-defined of FSLeyes <https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/customising.html#layouts>`_.

.. _Dicomifier: https://dicomifier.readthedocs.io
.. _FSLeyes: https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/
.. _spire: https://github.com/lamyj/spire
.. _WASSR: https://doi.org/10.1002/mrm.21873
