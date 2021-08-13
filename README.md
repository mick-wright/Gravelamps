# Gravelamps

## What is Gravelamps?

Gravelamps is a python package built upon the [`bilby`](https://git.ligo.org/lscsoft/bilby) framework and is designed to perform template based analysis of both simulated and real gravitationally lensed gravitational wave signals to determine the mass density profile of the lensing object. In so doing, it will also give estimates of the lens and source parameters for each of the mass density profiles that is tested. 

It is able to do this in both the wave optics and the geometric optics regimes - corresponding to lower and higher mass lenses as well as being able to cross between regimes in a single run allowing for the entire mass spectrum to be covered without break. The wave optics calculations represent particularly complex calculations which to be performed effectively are done using the C arbitrary precision library [`arb`](https://github.com/fredrik-johansson/arb/) with the overall lensing generation performed in a variety of C++ subprograms inside of Gravelamps.

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

Once installed, the process of performing an analysis run is very easy. Simply modify the INI file for the flavour of Gravelamps you wish to use (local for running on a local machine, pipe for running on a clusterised machine using the HTCondor Scheduler), and call the appropriate inference function giving it the INI as the only argument, these being

+ `gravelamps_local_inference` for local machine runs
+ `gravelamps_pipe_inference` for clusterised machine runs using the HTCondor Scheduler

In the case of `gravelamps_local_inference`, this will immediately start the run using the local resources available. In the case of `gravelamps_pipe_inference` this will generate a series of HTCondor submit files and an overarching DAG that can be submitted to the scheduler. In the latter case, the output of the program will detail the submission command.

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

#### Python Dependencies

In addition to a working `python 3` environment, the following are required for the overall running of the python package:
+ `numpy`
+ `scipy`
+ `astropy`
+ `bilby`

The `lalsuite` optional dependency of `bilby` is also required.

# Citation

The source code of Gravelamps is provided under the MIT License. If used for scientific purposes, please cite:

+ Paper forth coming

In addition, due to the dependency of Gravelamps upon in particular `bilby` and `arb`, but all of the mentioned components, please follow their citation practices. 
