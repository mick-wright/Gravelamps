# Gravelamps

Gravelamps is a model selection pipeline designed to determine the lens mass density profile from arbitrary lensed Gravitational Wave signals.

\| [Documentation](https://mick.wright.docs.ligo.org/Gravelamps) 
\| [Releases](https://git.ligo.org/mick.wright/Gravelamps/-/releases)
\| [Issue Tracker](https://git.ligo.org/mick.wright/Gravelamps/-/issues)
\| [Conda Page](https://anaconda.org/conda-forge/gravelamps)
\| 

![conda_version](https://anaconda.org/conda-forge/gravelamps/badges/version.svg) 
![pipeline](https://git.ligo.org/mick.wright/Gravelamps/badges/o4-development/pipeline.svg)
[![DOI](https://zenodo.org/badge/328470267.svg)](https://zenodo.org/badge/latestdoi/328470267)

Gravelamps is a python package built upon the [`bilby`](https://git.ligo.org/lscsoft/bilby) parameter estimation framework designed to perform analysis of gravitationally lensed gravitational wave signals to determine the most probable mass density profile of the lensing object. In so doing, it also aims to give estimates of the lens and source parameters for each of the mass density profiles tested.

Gravelamps is able to calculate the amplification factor for a variety of models in both the wave optics and geometric optics regimes. It allows the user to switch between these environments to allow a great deal of flexibility in the mass spectrum of the lensing object. In the wave optics regime, owing to the great difficulty of the calculations, Gravelamps makes use of the C arbitrary precision library [`arb`](https//arblib.org) and grants the choice of precision to the user to increase its suitability for the user's needs.

# Installation

Gravelamps is available via the conda package manager, and can be installed by running

	$ conda install -c conda-forge gravelamps

# Running Gravelamps

TO DO: create an example

# Citing Gravelamps

The source code of Gravelamps is provided under the MIT License. Citation of the code may be done using the DOI at the top

If Gravelamps has been useful for the production of scientific work, a citation for the methodological paper contained within Gravelamps is provided below:

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
	}

# Contributing to Gravelamps

The developers of Gravelamps welcome additional contribution. We particularly encourage the addition of lens models should you have need to write one. 

TO DO: Write a Contribution Guide
