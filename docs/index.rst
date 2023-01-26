===================
Gravelamps Overview
===================

        **Note**

        Gravelamps' Documentation is a living record and as such, you may find things that have broken over time or that have not been updated. We welcome contributions to this documentation, via the contribution guidelines outlined at the `Gravelamps Repo <https://git.ligo.org/mick.wright/Gravelamps>`_. 

Gravelamps
==========

Gravelamps is a python package built upon the `bilby <https://git.ligo.org/lscsoft/bilby>`_ framework designed to perform template based analysis of both simulated and real gravitationally lensed gravitational wave signals to determine the mass density profile of the lensing object. In so doing, it will also give estimates of the lens and source parameters for each of the mass density profiles that is tested. 

It is able to do this in regimes considering either only full wave optics, or considering both full wave optics as well as the simpler geometric optics case at a specified threshold allowing a great deal of flexibility in the mass spectrum of the lensing object. The particularly complex calculations required to compute the wave optics case are doing using the C arbitrary precision library `arb <https://arblib.org>`_ with the C++ libraries that have been written to implement these computations contained within ``gravelamps.model``.

Citation Guidelines
===================

.. image:: https://zenodo.org/badge/328470267.svg
   :target: https://zenodo.org/badge/latestdoi/328470267

The source code of Gravelamps is provided under the MIT License. If the software has been helpful to you, citation of the codemay be done using the DOI above. If used for scientific purposes, a citation is provided for the paper describing the methodology is provided below::

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
	journal = {The Astrophysical Journal},
	abstract = {We present the package Gravelamps, which is designed to analyze gravitationally lensed gravitational wave signals in order to constrain the mass density profile of the lensing object. Gravelamps does this via parameter estimation using the framework of bilby, which enables estimation of both the lens and the source parameters. The package can be used to study both microlensing and macrolensing cases, where the lensing mass distribution is described by a point-mass and extended-mass density profile, respectively. It allows the user to easily and freely switch between a full wave optics and approximate geometric optics description. The performance of Gravelamps is demonstrated via simulated analysis of both microlensing and macrolensing events, illustrating its capability for both parameter estimation and model selection in the wave optics and hybrid environments. To further demonstrate the utility of the package, the real gravitational-wave event GW170809 was analyzed using Gravelamps; this event was found to yield no strong evidence supporting the lensing hypothesis, consistent with previously published results.}
	}

In addition, due to the dependency of Gravelamps upon in particular, ``bilby``, ``bilby_pipe``, and ``arb``, but all of the mentioned components, please follow their citation practices.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   generate-lens
   inference
   ini

.. toctree::
   :maxdepth: 1
   :caption: Indices

   API Reference <autoapi/gravelamps/index>

