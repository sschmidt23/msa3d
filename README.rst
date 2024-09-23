MSA3D
=====


About
-----

This software was developed for data reduction and cube design for the JWST slit-stepping survey GO-2136.
It is designed to be applicable for any future JWST slit-stepping surveys employing a similar observing strategy.

The software is uses version 1.14.0 ``jwst`` STScI pipeline to process the data
via a modified set of arguments and keywords and original software dedicated to
cube design for a slit-stepping strategy with NIRSpec MSA --- an unsupported
processing mode in the standard STScI pipeline.  See  `Barisic et al. 2024
<https://ui.adsabs.harvard.edu/abs/2024arXiv240808350B/abstract>`__ . for
technical details and a case study analysis of an example target.


Installation
------------

To install MSA3D, clone the git repository:

.. code-block:: console

    git clone https://github.com/barisiciv/msa3d.git

which pulls the files into a directory called 'msa3d'.  To install, first create
a fresh python environment using conda or venv, then run:

.. code-block:: console

    pip install -e .

Disk space
----------

Total disk space required for full reduction (excluding STScI/Spec1Pipeline) is ~80GB, of which approximately:

    - 11GB : \*rate.fits files (available for download on MAST)

    - 50GB : products of custom JWST/STScI Spec2Pipeline + Spec3Pipeline reduction (2D spectra)

    - 12GB : products of cube design (data cubes and related products)


Data access
-----------

The ``msa3d`` data reduction starts with slope images, which can be found on
`MAST Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`__.
GO-2136 program is publicly available. Search for the data set on the portal
using the Proposal ID: 2136 and download all the *_rate.fits files.

After downloading the \*_rate.fits files, make sure all \*_rate.fits files are
in the **same folder**.


Running the software
---------------------

Temporary : **see notebook example**

