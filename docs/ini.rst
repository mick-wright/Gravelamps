=================
INI Configuration
=================

Both the Lens Generation and Parameter Estimation parts of Gravelamps are largely configured via an INI file. This INI file is compatible with both programs with the Lens Generation only INI being a subset of the options for a full parameter estimation INI. Below each option will be explained in some detail with suggested defaults. Settings used only by the inference program will be noted with a *. 

        **A Word of Caution**

        Where used, the section headers should be kept in the INI file. The Gravelamps INI is parsed using the python builtin ``configparser`` and consequently these section headers are used in the parsing.

Output Settings
===============

Example Rendering
-----------------

.. code-block:: ini
   
   [output_settings]
   outdir = example-output-directory
   label = example-run
   logging_level = INFO

Options
-------

* ``outdir`` : Name of the output directory into which all output subdirectories are placed. This option will default to the folder from which the program is run.
* ``label``: Name of the run as used in ``bilby`` based functions and applied to results from these functions.*
* ``logging_level``: Level of verbosity that will be printed by the Gravelamps logger. Possible options are ``INFO``, ``WARNING``, and ``CRITICAL``. This setting will be overridden in the program is run with ``-v`` or ``--verbose`` options. The default is ``INFO``.

Run Settings
============

Example Rendering
-----------------

.. code-block:: ini

   [run_settings]
   local = False
   injection = True
   submit = False
   millilensing_kmax = 2

Options
-------

* ``local``: Flag to determine whether the run will be performed locally or be built for the HTCondor scheduler. Defaults to False.
* ``injection``: Flag to determine whether the run contains an injection setup. Defaults to False.
* ``submit``: Flag to determine whether when built for the HTCondor scheduler, the run should be automatically submitted. Defaults to False.
* ``millilensing_kmax``: The maximum number of millilensing images that will be considered in an inference run using the millilensing framework.*

Condor Settings
===============

Example Rendering
-----------------

.. code-block:: ini

   [condor_settings]
   request_cpus = 16
   request_memory = 8 GB
   request_disk = 5 GB
   accounting_group = ligo.dev.o4.cbc.lensing.multi

Options
-------

* ``request_cpus``: How many CPU cores to request from the scheduler.
* ``request_memory``: How much RAM to request from the scheduler.
* ``request_disk``: How much storage space to request from the scheduler.
* ``accounting_group``: Accounting tag to use for the scheduler.

Additional Options
------------------

The above options are the necessary settings for a run to be built for the HTCondor scheduler, however, HTCondor has a variety of potential options, all of which may be set from within this section using the syntax provided by the `HTCondor User Manual<https://htcondor.readthedocs.io/en/latest/>_`  

Lens Generation Settings
=========================

Example Rendering
-----------------

.. code-block:: ini

   [injection_lens_generation_settings]
   lensing_module = gravelamps.lensing.interpolator
   interpolator_model = point

   minimum_dimensionless_frequency = 0.01
   maximum_dimensionless_frequency = 1000
   length_dimensionless_frequency = 10000

   minimum_source_position = 0.1
   maximum_source_position = 3.0
   length_source_position = 20

   arithmetic_precision = 1024
   geometric_optics_frequency = 1001

   [analysis_lens_generation_settings]
   lensing_module = gravelamps.lensing.interpolator
   interpolator_model = sis

   dimensionless_frequency_file = /path/to/dimensionless-frequency.dat
   source_position_file = /path/to/source-position.dat
   amplification_factor_real_file = /path/to/amplification-factor-real.dat
   amplification_factor_imag_file = /path/to/amplification-factor-imag.dat

Options
-------

Settings below may go in either section and refer to the injection or analysis waveform respectively.

* ``lensing_module``: Python path to a module capable of generating amplification factor data. Within Gravelamps these are contained within ``gravelamps.lensing``.
* ``interpolator_model``: When using the ``gravelamps.lensing.interpolator`` module, this instructs the module which model of lens is being used.

The grid parameters of dimensionless frequency and source position are handled both in the same manner as:

* ``minimum_{parameter}``: Minimum value of the grid parameter
* ``maximum_{parameter}``: Maximum value of the grid parameter
* ``length_{parameter}``: Number of grid points to generate between minimum and maximum

Next are the parameters required when using the ``gravelamps.lensing.interpolator`` module which allows access to both wave and geometric optics functions:

* ``arithmetic_precision``: Bit precision to use for wave optics calculations
* ``geometric_optics_frequency``: Dimensionless frequency value at which to switch from wave optics to geometric optics. May be set at 0 to always use geometric optics or greater than the maximum dimensionless frequency to always use wave optics.
* ``sis_upper_summation_limit``: Needed when ``ìnterpolator_model`` is ``sis``. Determines the maximum term of the summation calculated in the amplification factor calculations.
* ``nfw_upper_integration_limit``: Needed when ``interpolator_model`` is ``nfw``. Determines the finite upper bound of the integration calculated in the amplification factor calculations. 
* ``nfw_scaling_constant``: Needed when ``ìnterpolator_model`` is ``nfw`` or when ``lensing_module`` is ``gravelamps.lensing.nfw``. Determines the characteristic scale of the NFW profile.

The final set of options may be used with ``gravelamps.lensing.interpolator`` to signify that the data is already computed:

* ``dimensionless_frequency_file``: Path to file containing dimensionless frequency grid points
* ``source_position_file``: Path to file containing source position grid points
* ``amplification_factor_real_file``: Path to file containing grid of the real component of the amplification factor.
* ``amplification_factor_imag_file``: Path to file containing grid of the imaginary component of the amplification factor.

The amplification factor files used in this section should be arranged as matrices of source position by dimensionless frequency. If they are not square, incorrectly oriented files can be transposed, however, this cannot be done for square matrices as there is no way to verify which orientation is correct.

Inference Settings*
===================

Example Rendering
------------------

.. code-block:: ini

   [inference_settings]
   detectors = H1, L1, V1
   duration = 4.0
   sampling_frequency = 1024
   trigger_time = 0.0
   sampler = dynesty
   prior-file = /path/to/file
   waveform-generator-class = gravelamps.lensing.waveform_generator.LensedWaveformGenerator
   frequency-domain-source-model = lal_binary_black_hole
   sampler-kwargs = {'nlive':1000, 'naccept':60, 'check_point_plot':True, 'check_point_delta_t':1800, 'print_method':'interval-60', 'sample':'acceptance-walk'}

   waveform_approximant = IMRPhenomXPHM
   reference_frequency = 20
   minimum_frequency = 20
   maximum_frequency = 448

Options
-------

* ``detectors``: List of the interferometers in which the signal has been detected.
* ``duration``: Duration of the signal in seconds
* ``sampling_frequency``: Sampling rate of the signal in Hz
* ``trigger_time``: GPS time at which the signal trigger was identified
* ``sampler``: Nested sampler to be used. Any sampler accepted by ``bilby`` may be used, ``dynesty`` is recommended.
* ``prior-file``: File containing the priors for the run.
* ``waveform-generator-class``: Python path to the waveform generator to be used. Will default to ``gravelamps.lensing.waveform_generator.LensedWaveformGenerator``, though any derivative class of ``bilby.gw.waveform_generator.WaveformGenerator`` may be used.
* ``frequency-domain-source-model``: Python path to the frequency domain source model. Defaults to ``lal_binary_black_hole``.
* ``sampler-kwargs``: Dictionary containing options that will be passed directly to the sampler.
* ``waveform_approximant``: Numerical relativity approximant that will be used. 
* ``reference_frequency``: Reference frequency for the waveform in Hz.
* ``minimum_frequency``: Minimum frequency of the waveform to be generated in Hz.
* ``maximum_frequency``: Maximum frequency of the waveform to be generated in Hz.

Bilby Pipe Settings*
====================

Example Rendering
-----------------

.. code-block:: ini

   [bilby_pipe_additional_settings]
   n-simulation = 1
   n-parallel = 2
   generation-seed = 1234

Options
-------

Settings here will be passed to `bilby_pipe<https://lscsoft.docs.ligo.org/bilby_pipe/master/index.html>_` and should follow the appropriate syntax. This is also where settings such as ``psd-dict`` that govern GW data fetching should be placed, again following the structure of ``bilby_pipe``.

Injection Parameters*
=====================

Example Rendering
-----------------

.. code-block:: ini

   [injection_parameters]
   mass_1 = 36.0
   mass_2 = 29.0
   a_1 = 0.4
   a_2 = 0.3
   tilt_1 = 0.5
   tilt_2 = 1.0
   phi_12 = 1.7
   phi_jl = 0.3
   luminosity_distance = 410
   theta_jn = 0.4
   phase = 1.3
   ra = 1.375
   dec = 1.12108
   psi = 2.659
   lens_mass = 1000
   source_position = 0.3
   geocent_time = 1125259642.413
   lens_fractional_distance = 0.5

Options
-------

These largely follow the LAL parameters needed for the source model. 

The Gravelamps specific parameters for the microlensing models are as follows:

* ``lens_mass``: Source frame (non-redshifted) mass of the lensing object. 
* ``source_position``: Dimensionless displacement from the optical axis.
* ``lens_fractional_distance``: Distance to the lensing object as a fraction of the luminosity distance to the source.

For the model-agnostic millilensing model, the number of parameters will vary. Those parameters that vary will be written as ``{parameter}x`` where ``x`` indicates that one parameter will be needed for each image beyond the first.

* ``k``: Number of images in the millilensing signal.
* ``dlx``: Effective luminosity distance of the image signal.
* ``dtx``: Time delay between signal relative to previous image.
* ``n0``: Morse phase of the first signal.
* ``nx``: Morse phase of the image.
