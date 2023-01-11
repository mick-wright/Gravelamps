# Gravelamps

## What is Gravelamps?

Gravelamps is a python package built upon the [`bilby`](https://git.ligo.org/lscsoft/bilby) framework designed to perform template based analysis of both simulated and real gravitationally lensed gravitational wave signals to determine the mass density profile of the lensing object. In so doing, it will also give estimates of the lens and source parameters for each of the mass density profiles that is tested. 

It is able to do this in both the wave optics only regime as well as in a hybrid environment crossing into the geometric optics regime at a specified threshold allowing a great deal of flexibility in the mass spectrum of the lensing object. The particularly complex calculations required to compute the wave optics case are done using the C arbitrary precision library [`arb`](https://github.com/fredrik-johansson/arb/) with the C++ libraries that have been written to implement these computations are contained within `gravelamps.model`. 

## How Does It Work?

The anatomy of a Gravelamps analysis run is as follows:

1. Get/Generate a set of dimensionless frequencies and source positions over which to consider
2. Construct the value of the amplification factor for a given lens profile over these ranges
3. Generate an amplification factor interpolator over the area of interest
4. Construct a lensed waveform generator based upon this amplification factor
5. Get/Generate the data to be estimated
6. Perform a parameter estimation run using bilby for the profile

This process will yield estimates of the lens and source parameters as well as the evidence for that lens profile. By performing multiple runs with multiple profiles, these evidences can be compared to give the preferred overall lensing model. 

## How Do I Perform A Run?

Once installed, the process of performining an alysis run is very easy. Simply modify the example INI file given to the specifications that are desired for your analysis run and call the inference program as follows:

`gravelamps_inference yourfile.ini`

This will then either run the jobs directly or generate a DAG which you can either submit yourself or have Gravelamps submit for you depending upon settings within the INI. Additional help may be gathered by running:

`gravelamps_inference -h`

Should you only desire to have Gravelamps perform the generation of data for a lensing interpolator, this may be done by running

`gravelamps_generate_lens yourfile.ini`

## Installation Instructions

The process of installation has designed to be as simple as possible. Assuming all necessary dependencies are met, simply follow this procedure:

1. Clone the repository
2. In a terminal, navigate to the repository directory
3. `pip install .` 

This can be done either in default python environment or in a conda environment without issue. The C++ subprograms will by default be placed into $HOME/bin - which the user should have extant already and on the $PATH environment variable so that the calls can be made without issue. The behaviour of the C++ subprogram installation can be modified via the provided makefile in the lensing subfolder - which will be called when the pip installation is made. 

## Dependencies

#### C++ Dependencies

The following are required for the compilation and running of the C++ subprograms for lens generation:
+ `arb`
+ `openMP`
+ `g++`
+  `GNU make`
+  `Boost`

#### Python Dependencies

In addition to a working `python 3` environment, the following are required for the overall running of the python package:
+ `numpy`
+ `scipy`
+ `astropy`
+ `bilby`
+ `htcondor`
+ `bilby_pipe`

The `lalsuite` optional dependency of `bilby` is also required.

# Citation

[![DOI](https://zenodo.org/badge/328470267.svg)](https://zenodo.org/badge/latestdoi/328470267)

The source code of Gravelamps is provided under the MIT License. If the software has been helpful to you, citation of the code may be done using the DOI above. If used for scientific purposes, a citation is provided for the paper describing the methodology is provided below:
	
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

In addition, due to the dependency of Gravelamps upon in particular `bilby`, `bilby_pipe`, and `arb`, but all of the mentioned components, please follow their citation practices. 

# Contribution Guidelines

The developers welcome additional contribution - one of the major design intentions for Gravelamps is to make an easily extensible platform - adding additional lens models and scenarios as and when users need them, and we encourage users to share these models when they are made. To make sure everything goes smoothly, if you wish to contribute to Gravelamps, please follow the following procedures:

1. Open an issue detailing the change that you model that you intend to implement, ideally providing a timescale, this will minimise the amount of work duplication that occurs. 
2. Fork the repository to allow you to work on your own branch
3. Create the model as well as any additional work that needs doing to implement. 
4. Code should be linted, that coming from the authors is linted using pylint for the python sections and using cpplint for the C++ sections
5. Once complete and ready to be added to Gravelamps fully, issue a merge request citing your issue.
