===================
Gravelamps Overview
===================

.. image:: https://img.shields.io/badge/Anaconda.org-2.3.1-blue.svg?style=flat-round
   :target: https://anaconda.org/conda-forge/gravelamps
.. image:: https://git.ligo.org/mick.wright/Gravelamps/badges/o4-development/pipeline.svg
   :target: https://git.ligo.org/mick.wright/Gravelamps

Gravelamps
==========

Gravelamps is a python package built upon the `bilby <https://git.ligo.org/lscsoft/bilby>`_ parameter estimation framework designed to perform analysis of gravitationally lensed gravitational wave signals to determine the most probable mass density profile of the lensing object. In so doing, it also aims to give estimates of the lens and source parameters for each of the mass density profiles tested.

Gravelamps is able to calculate the amplification factor for a variety of models in both the wave optics and geometric optics regimes. It allows the user to switch between these environments to allow a great deal of flexibility in the mass spectrum of the lensing object. In the wave optics regime, owing to the great difficulty of these calculations, Gravelamps makes use of the C arbitrary precision library `arb <https://arblib.org>`_ and grants the choice of precision to the user to increase its suitability for the user's needs.

Citing Gravelamps
=================

.. image:: https://zenodo.org/badge/328470267.svg
   :target: https://zenodo.org/badge/latestdoi/328470267

Gravelamps is freely provided under the MIT License. If it is helpful to you, please consider citing the code which may be done using the above DOI. If Gravelamps has been helpful specifically to work heading to publication, a citation of the methodological paper presenting Gravelamps is provided below:

.. code-block:: bibtex

	@article{Wright_2022,
	         doi = {10.3847/1538-4357/ac7ec2},
	         url = {https://doi.org/10.3847/1538-4357/ac7ec2},
	         year = 2022,
	         month = {aug},
	         publisher = {American Astronomical Society},
	         volume = {935},
	         number = {2},
	         pages = {68},
	         author = {Mick Wright and Martin Hendry},
	         title = {Gravelamps: Gravitational Wave Lensing Mass Profile Model Selection},
	         journal = {The Astrophysical Journal}
	        }

In addition, Gravelamps is very much built on top of existing frameworks, in particular ``bilby``, ``bilby_pipe``, and ``arb``. Please also consider following the citation guidelines of these and all of Gravelamps' dependencies.

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   installation
   generate-lens
   inference
   ini

.. toctree::
   :maxdepth: 1
   :caption: Code Index

   API Reference <autoapi/gravelamps/index>

