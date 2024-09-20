# msa3d



About
------------

This software was developed for data reduction and cube design for the JWST slit-stepping survey GO-2136.
It is designed to be applicable for any future JWST slit-stepping surveys employing a similar observing strategy.

The software is uses a combination of a modified version of the standard ``jwst`` STScI pipeline, and the original software dedicated to cube design for a slit-stepping strategy with NIRSpec MSA -- an unsupported processing mode in the standard STScI pipeline.
See  `Barisic et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024arXiv240808350B/abstract>`_ . for technical details and a case study analysis of an example target.


Installation
------------

To install MSA3D, clone the git repository:

.. code-block:: console

    git clone https://github.com/barisiciv/msa3d.git

Which pulls the files into a directory called 'msa3d'

In order to execute MSA3D code, create a new conda environment using the specifications in the ``environment.yml`` file.
This software uses 1.14.0 version of the ``jwst`` package (other dependencies specified in the ``environment.yml`` file).


Disk space
------------

Total disk space required for full reduction (excluding STScI/Spec1Pipeline) is ~80GB, of which approximately:

    - 11GB : *rate.fits files (available on MAST)

    - 50GB : products of custom JWST/STScI Spec2Pipeline + Spec3Pipeline reduction (2D spectra)

    - 12GB : cubes and related data products


Data access
------------

The ``msa3d`` data reduction starts with slope images, which can be found on `MAST Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`__. GO-2136 program is publicly available. Search for the data set on the portal using the Proposal ID: 2136 and download all the *_rate.fits files.

After downloading, make sure all *_rate.fits files are in the **same folder**.


Running the software
------------

