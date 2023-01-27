========================
Gravelamps INI Structure
========================

The following documents the INI that is used as the input to the Gravelamps programs &mdash; detailing what each feature does as well as their usage within Gravelamps.

         **A Word of Caution**

         The sections as specified below within the Gravelamps INI should be maintained where the settings are applicable, the Gravelamps INI is parsed using the builtin ``configparser`` and consequently the keys are referenced by section. 

Output Settings
===============

        Rendered as ``[output_settings]``.

Settings within this section, as might be inferred from the name, concern the location and verbosity of output from the Gravelamps programs. Settings are applicable to all Gravelamps programs.

* ``outdir`` : Name of the folder into which the output subdirectories are placed. Defaults to the folder from which the program is run. 
* ``label``: *Applicable to programs that use ``bilby`` functions only* The naming label associated with the run within ``bilby``.
* ``logging_level``: Level of verbosity that will be printed by the Gravelamps logger. Settings may be ``INFO``, ``WARNING``, ``CRITICAL``. Setting may be overridden by running the program with the ``--verbose``/``-v`` option. Will default to ``INFO`` level. 

Run Settings
============

	 Rendered as ``[run_settings]``

Settings within this section are INI setting versions of the flags that may be set when running the programs, and may for each option be also set by running the program including ``--option_name``.

* ``local``: If set to True, all parts of the run will be performed locally as opposed to using the ``HTCondor`` scheduler. Defaults to False.
* ``injection``: If set to True, the run contains an injection setup. Defaults to False. 
* ``submit``: If set to True, the Condor submission files will be directly submitted to the scheduler. 

Condor Settings
===============

	 Rendered as ``[condor_settings]``

Settings placed here are options that will be placed within any generated ``HTCondor`` submit files created by the programs. 

This section is optional when the run is set to local. 

For a list of possible options, consult the ``HTCondor`` documentation. A warning will be given by the programs should they attempt to be ran without the ``accounting_group`` value being set. 

Lens Generation Settings
========================

	 Rendered as ``[analysis_lens_generation_settings]`` and ``[injection_lens_generation_settings]`` for the analysis and injection waveforms respectively. 

These settings may be added to both sections. The section for the analysis lensing generation is always necessary whilst the section for the injection lensing generation is necessary only where an injection is needed &emdash; noted either by the setting ``injection = True`` in the ``[run_settings]`` section, or by the use of the argument flag ``--injection``. 

* ``lensing_module``: Python path to where the lensing model that is being selected is, this may be any available python path. However, those implemented within Gravelamps are placed within ``gravelamps.lensing``. Available choices are ``interpolator``, ``point``, ``sis``, and ``nfw``. 
* ``interpolator_model``: If using ``gravelamps.lensing.interpolator`` above, this specifies if the interpolator is adhering to one of the implemented Gravelamps models. This allows for generation of interpolator data. If it is left blank, no new data can be constructed though inference may be performed.

For each of ``dimensionless_frequency`` and ``source_position``, the following settings determine how to construct the grid of these parameters should new interpolator data be being generated.
* ``minimum_{quantity_name}``: Minimum value of the quantity
* ``maximum_{qunatity_name}``: Maximum value of the quantity
* ``length_{quantity_name}``: How many points to generate between the minimum and maximum

The following settings determine how the calculations will be performed when new data is being generated. Settings below that are not preceded by a lens model are necessary for all lens models. Those preceded by a lens model are necessary only for that particular model.

* ``arithmetic_precision``: Bit precision to use for wave optics calculations.
* ``geometric_optics_frequency``: Dimensionless frequnecy at which to switch from wave optics calculations to geometric optics calculations. Lower values will trade speed of calculation for accuracy of calculation. 

* ``sis_upper_summation_limit``: Upper limit of the summation used to calculate the SIS wave optics amplification factor.
* ``nfw_upper_integration_limit``: Upper limit of the integration used to calculate the NFW wave optics amplification factor.
* ``nfw_scaling_constant``: Characteristic scaling of the NFW profile. Used in both wave and geometric optics calculations.

The following settings are optional, and when used will bypass the generation of new data in favour of using the data from these files. If the files are not a complete set, the other options will be used to complete them.

* ``dimensionless_frequency_file``: File containing dimensionless frequency data. 
* ``source_position_file``: File containing source position data.
* ``amplification_factor_real_file``: File containing real component of the amplification factor.
* ``amplification_factor_imag_file``: File containing imaginary component of the amplification factor.

Amplification factor files used for this section should be arranged as matrices of source position by dimensionless frequency, or vice versa. If not square these will be corrected if arranged incorrectly, however, they must be corrected oriented if square. 

Inference Settings
==================

	 Rendered as ``[inference_settings]``

Settings within this section concern how the waveforms should be generated and settings that are necessary for the construction/fetching of the data. 

* ``detectors``: List of Interferometers in which the signal has been detected.
* ``duration``: Duration of the signal in seconds.
* ``sampling_frequency``: Sampling frequency for the signal.
* ``trigger_time``: GPS Time at which the trigger for the signal occurred.
* ``sampler``: Nested sampler to be used for the inference. Any sampler that may be used by ``bilby`` is accepted.
* ``prior-file``: File containing the priors for the run. See ``bilby`` documentation for more information on constructing these files.
* ``waveform-generator-class``: Python path containing the waveform generator to be used. Gravelamps' own generator is ``gravelamps.lensing.waveform_generator.LensedWaveformGenerator`` and is the most supported, but any child class of ``bilby.gw.waveform_generator.WaveformGenerator`` will be accepted.
* ``frequency-domain-source-model``: Source model to be used. Defaults to LAL's Binary Black Hole.
* ``sampler_kwargs``: Dictionary of arguments that can be passed directly to the sampler.
* ``waveform_approximant``: Approximant for the waveform to be used.
* ``reference_frequency``: Reference frequency for the waveform
* ``minimum_frequency``: Minimum frequency of the waveform
* ``maximum_frequency``: Maximum frequency of the waveform

Bilby Pipe Additional Settings
==============================

	 Rendered as ``[bilby_pipe_additional_settings]``

Any settings placed here will be passed to ``bilby_pipe`` functions as part of the inference process. These must be given as they would in the ``bilby_pipe`` INI and are documented within the ``bilby_pipe`` documentation. One should place settings such as ``n-simulation`` for injections or the ``channel-dict`` and ``psd-dict`` for actual GW data as adhering to the ``bilby_pipe`` standards. 

Injection Parameters
====================

	 Rendered as ``[injection_parameters]``

Parameter specifications for the injected signal. For full documentation on this see ``bilby`` documentation.
