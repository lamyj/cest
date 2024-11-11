# CEST Toolbox

## Installation

The CEST toolbox is based on Python packages ([doit](https://pydoit.org), [nibabel](https://nipy.org/nibabel), [numpy](https://numpy.org), [scipy](https://scipy.org), [spire](https://github.com/lamyj/spire), [pyyaml](https://pyyaml.org)) and additional software ([ANTs](https://github.com/ANTsX/ANTs), [Dicomifier](https://dicomifier.readthedocs.io), [FSLeyes](https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/)).

The recommended way is to install the CEST toolbox through an Anaconda-compatible environment, e.g. [Miniforge](https://github.com/conda-forge/miniforge). With this, the dependencies can be installed with  `conda install -c conda-forge ants dicomifier fsleyes nibabel numpy scipy spire-pipeline pyyaml`.

## Usage

There are two main ways to use this toolbox: through function calls, or through tasks in a [spire](https://github.com/lamyj/spire) pipeline. The former approach works well for experimenting and when processing a limited quantity of data while the latter approach is particularly suited for reproducible analysis on a larger database.

### Functions

Functions work on files: they read and write data from the disk, and thus have no return value. Images, including multi-volume images of Z-spectra, are expected to be stored in [well-known format](https://nipy.org/nibabel/api.html). Meta-data are expected in the JSON format used by [Dicomifier](https://dicomifier.readthedocs.io).

For example, a shift map using the [WASSR](https://doi.org/10.1002/mrm.21873) method can be computed using the following code sample. It will load data from a Z-spectrum map and its associated meta-data, and store the result in `deta_ppm.nii.gz`.

```python
import cest

cest.wassr("z_spectrum.nii.gz", "z_spectrum.json", "delta_ppm.nii.gz")
```

The documentation of individual functions is described in the [API documentation](docs/api/functions.rst).

### Tasks

By linking tasks in a [spire](https://github.com/lamyj/spire) pipeline, the depencies between processing steps are handled automatically, improving the reproducibility of your pipeline. The example below assumes that individual exams are stored as subdirectories of `/storage/my_study`, and, for each exam, creates a shift map in PPM, shifts the Z-spectra and resamples it in the [-5, +5] PPM range.

```python

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
```

The documentation of tasks is available in the [API documentation](docs/api/tasks.rst).

### Visualization

[FSLeyes](https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/) can load NIfTI images and display the Z-spectrum of a single voxel. The provided [layout file](docs/_static/cest.txt) can be [added to the pre-defined of FSLeyes](https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/customising.html#layouts).
