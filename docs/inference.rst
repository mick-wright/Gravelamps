===========================
Lensed Parameter Estimation
===========================

One of the major purposes of Gravelamps is to be a framework for lens model selection of the implemented models. One of the primary ways of doing this is to perform parameter estimation of arbitray lensed gravitational wave signals. It does this by integrating with the parameter estimation framework ``bilby`` to gain access to its wide array of gravitational wave signal waveforms and nested samplers. 

The program which performs a parameter estimation run using the Gravelamps framework is ``gravelamps_inference`` and may typically be run as follows::

        $ gravelamps_inference -o --option /path/to/ini-file.ini

with a list of the possible options and an example INI configuration with explanation given below. Note that the parameter estimation program is fully complete, meaning that it is capable of running ``gravelamps_generate_lens`` by itself to simplify the processs of running for users. 

Inference Options
=================

-h, --help              Shows the help message and exits
-i, --injection         Run contains an injection and so injection settings should be considered as well
-l, --local             Run should be performed locally rather than built for the HTCondor scheduler
-s, --submit            Run should be directly submitted to the HTCondor scheduler upon being built
-v, --verbose           Display the maximal level of information from the logger, overriding the INI settings

Example Inference
=================

The following example will perform parameter estimation for GW150914-like gravitational wave signal that has been lensed by an isolated point mass of 1000 solar masses at a dimensionless source position of 0.3 that has been placed half-way between the source and the observer. 

To run the example, copy the example configuration below to a file called ``example.ini`` and the example prior configuration to a file called ``example.prior`` and run the inference program as follows::

        $ gravelamps_inference example.ini

This will generate a folder called ``inference-example`` which will contain all of the submission and result data. To submit the run to the HTCondor scheduler so as to perform all of the necessary analysis, run::

        $ condor_submit_dag inference-example/submit/gravelamps_inference.dag

Once complete, the results will be contained within the folders. Any generated amplification factor data will be stored in the ``inference-example/data`` folder and the final results of the inference investigation will be stored in ``inference-example/final_result`` with the component chain results in ``inference-example/result`` as they would be with any ``bilby`` run.

Example INI
-----------

.. code-block:: ini

        [output_settings]
        # Sets the output directory to 'inference-example'
        outdir = inference-example
        # Sets the label bilby gives to the run to lensed-ex
        label = lensed-ex
        # Sets default logging level
        logging_level = INFO

        [run_settings]
        # Run should be built for the HTCondor scheduler
        local = False
        # Run should not be submitted immediately upon being built
        submit = False
        # Run contains an injected signal
        injection = True

        [condor_settings]
        # Run requests 16 cores from the scheduler
        request_cpus = 16
        # Run requests 12 GB of RAM from the scheduler
        request_memory = 12 GB
        # Run requests 5 GB of storage space from the scheduler
        request_disk = 5 GB
        # Accounting tag for the scheduler, this is suitable for use on the LDG clusters
        accounting_group = ligo.dev.o4.cbc.lensing.multi

        [analysis_lens_generation_settings]
        # These settings apply to the analysis waveform which will be a point mass microlensing waveform
        # Set the lensing module to interpolator to use the wave optics calculations
        lensing_module = gravelamps.lensing.interpolator
        # Setting the lens model to the isolated point mass case
        interpolator_model = point

        # Dimensionless frequency grid settings
        minimum_dimensionless_frequency = 0.01
        maximum_dimensionless_frequency = 2000
        length_dimensionless_frequency = 100000

        # Source position grid settings
        minimum_dimensionless_frequency = 0.1
        maximum_source_position = 3.0
        length_source_position = 30

        # Arithmetic precision to use for wave optics calculations
        arithmetic_precision = 2048
        # Dimensionless frequency at which to switch to geometric optics
        geometric_optics_frequency = 1000

        [injection_lens_generation_settings]
        # These settings apply to the injection waveform whcih will be a point mass microlensing waveform
        # For this example these will be the same as the analysis settings, but can in principle be different

        lensing_module = gravelamps.lensing.interpolator
        interpolator_model = point

        minimum_dimensionless_frequency = 0.01
        maximum_dimensionless_frequency = 2000
        length_dimensionless_frequency = 100000

        minimum_source_position = 0.1
        maximum_source_position = 3.0
        length_source_position = 30

        arithmetic_precision = 2048
        geometric_optics_frequency = 1000

        [inference_settings]
        # Detectors which will be investigated for the signal
        detectors = H1, L1, V1
        # Duration of the signal in seconds
        duration = 4.0
        # Sampling frequency
        sampling_frequency = 1024
        # Trigger time for the GW event in GPS seconds
        trigger_time = 1125259642.413
        # Nested sampler to use
        sampler = dynesty
        # Location of the prior file, included below
        prior-file = example.prior
        # Waveform generator to use for the run
        waveform-generator = gravelamps.lensing.waveform_generator.LensedWaveformGenerator
        # Arguments to be passed to the sampler
        sampler-kwargs = {'nlive':1000, 'naccept':60, 'check_point_plot':True, 'check_point_delta_t':1800, 'print_method':'interval-60', 'sample':'acceptance-walk'}

        # NR-approximant to use
        waveform_approximant = IMRPhenomXPHM
        # Reference Frequency
        reference_frequency = 20
        # Waveform minimum frequency 
        minimum_frequency = 20
        # Waveform maximum frequency
        maximum_frequency = 1024

        [bilby_pipe_additional_settings]
        # Setting that one waveform should be injected
        n-simulation = 1
        # Allowing the waveform to enter parts of parameter space that would typically result in durations of signal longer than specified
        enforce-signal-duration = False
        # Noise seed into which the injection will be sent. Allows for reproducibility
        generation-seed = 1234
        # Number of parallel analysis chains to set off
        n-parallel = 2

        [injection_parameters]
        # True parameter values for the injection signal
        # Binary black hole source parameters
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
        geocent_time = 1125259642.413
        # Source frame mass of the lensing object in solar masses
        lens_mass = 1000
        # Dimensionless displacement from the optical axis
        source_position = 0.3
        # Fraction of the luminosity distance to place the lens at
        lens_fractional_distance = 0.5

Example Prior
-------------

.. code-block::

        chirp_mass = Uniform(name='chirp_mass', minimum=22, maximum=80, unit='$M_{\odot}$')
        mass_ratio = Uniform(name='mass_ratio', minimum=0.125, maximum=1)
        mass_1 = Constraint(name='mass_1', minimum=1.001398, maximum=1000)
        mass_2 = Constraint(name='mass_2', minimum=1.001398, maximum=1000)
        a_1 = Uniform(name='a_1', minimum=0, maximum=0.88)
        a_2 = Uniform(name='a_2', minimum=0, maximum=0.88)
        tilt_1 = Sine(name='tilt_1')
        tilt_2 = Sine(name='tilt_2')
        phi_12 = Uniform(name='phi_12', minimum=0, maximum=2 * np.pi, boundary='periodic')
        phi_jl = Uniform(name='phi_jl', minimum=0, maximum=2 * np.pi, boundary='periodic')
        luminosity_distance = PowerLaw(name='luminosity_distance', minimum=1e2, maximum=15000, alpha=2, unit='Mpc')
        dec = Cosine(name='dec')
        ra = Uniform(name='ra', minimum=0, maximum=2 * np.pi, boundary='periodic')
        theta_jn = Sine(name='theta_jn')
        psi = Uniform(name='psi', minimum=0, maximum=np.pi, boundary='periodic')
        phase = Uniform(name='phase', minimum=0, maximum=2 * np.pi, boundary='periodic')
        geocent_time = Uniform(minimum=1126259642.413-2, maximum=1126259642.413+2)
        source_position = PowerLaw(alpha=1, minimum=0.1, maximum=3.0, latex_label='$y$')
        lens_mass = Uniform(minimum=1, maximum=10000, latex_label='$M_{Lens}$', unit='$M_{\odot}$')
        lens_fractional_distance = 0.5
