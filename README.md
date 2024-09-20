# msa3d



About
------------

This software was developed for data reduction and cube design for the JWST slit-stepping survey GO-2136.
It is designed to be applicable for any future JWST slit-stepping surveys employing a similar observing strategy.

See  `Barisic et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024arXiv240808350B/abstract>`_ . for technical details and a case study analysis of an example target.


Installation
------------

To install MSA3D, clone the git repository:

.. code-block:: console

    git clone https://github.com/barisiciv/msa3d.git

Which will pull the files into a directory called 'msa3d'

In order to execute MSA3D code, create a new conda environment using the specifications in the ``environment.yml`` file.
This software uses ``jwst`` package version 1.14.0 (other dependencies specified in the ``environment.yml`` file).



Data access
------------

The ``msa3d`` data reduction starts with slope images, which can be found on `MAST Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`__. Search for the data set using 

