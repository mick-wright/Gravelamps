===============
Lens Generation
===============

Gravelamps is capable of constructing amplification factor data for a number of lens mass density profiles. It is able to construct this data for both the wave optics and geometric optics regimes and can cross between these two regimes at a user-specified point in order to increase utility to the user. As calculation of the amplification factor in the wave optics regime is computationally challenging particularly at higher dimensionless frequency, it is done using arbitrary precision with the user able to specify the precision that they need. 

Generation of amplification factor data is handled by the ``gravelamps_generate_lens`` program which is typically run as follows::

        $ gravelamps_generate_lens -o --option /path/to/ini-file.ini

with a list of the possible options and an example INI configuration with explanation given below.

Lens Generation Options
=======================

-h, --help              Shows the help message and exits
-i, --injection         Run is an injection and should use the injection instead of analysis lens settings
-l, --local             Run should be performed locally rather than built for the HTCondor scheduler
-s, --submit            Run should be directly submitted to the HTCondor scheduler upon being built
-v, --verbose           Display the maximal level of information from the logger, overrriding the INI settings

Example Lens Generation
=======================

The following example will generate the amplification factor for an isolated point mass purely using wave optics from a dimensionless frequency of 0.01 to a dimensionless frequency of 1000 for 20 source positions ranging from 0.1 to 3.0. This run will be built for the HTCondor scheduler. 

To run the example, copy the example configuration below to a file called ``example.ini`` and run the lens generation program as follows::

        $ gravelamps_generate_lens example.ini

This will generate a folder called ``lens-generation-example`` which will contain all of the submission and result data. To submit the run to the HTCondor scheduler so as to actually perform the lens generation, then run::

        $ condor_submit lens-generation-example/submit/generate_analysis_interpolator_data.sub

Once complete, the resultant amplification factor data can be found in the ``lens-generation-example/data`` folder. This will be in the form of four ``.txt`` files containing the dimensionless frequency and source position grid positions and the amplification factor's real and imaginary components.  

Example INI
-----------

.. code-block:: ini

        [output_settings]
        # Sets the output directory to 'lens-generation-example'
        outdir = lens-generation-example
        # Sets default logging level
        logging_level = INFO

        [run_settings]
        # Run should be built for the HTCondor scheduler
        local = False
        # Run should not be submitted immediately upon being built
        submit = False
        # Run is not an injection
        injection = False

        [condor_settings]
        # Run requests 16 cores from the scheduler
        request_cpus = 16
        # Run requests 8 GB of RAM from the scheduler
        request_memory = 8 GB
        # Run requests 2 GB of storage space from the scheduler
        request_disk = 2 GB
        # Accounting tag for the scheduler, this is suitable for use on the LDG clusters
        accounting_group = ligo.dev.o4.cbc.lensing.multi

        [analysis_lens_generation_settings]
        # Should use the interpolator module for access to the wave optics calculations
        lensing_module = gravelamps.lensing.interpolator
        # Setting the lens model to the isolated point mass case
        interpolator_model = point
        
        # Dimensionless frequency grid settings
        minimum_dimensionless_frequency = 0.01
        maximum_dimensionless_freuqency = 1000
        length_dimensionless_frequency = 10000

        # Source position grid settings
        minimum_source_position = 0.1
        maximum_source_position = 3.0
        length_source_position = 20

        # Arithmetic precision to use for wave optics calculations
        arithmetic_precision = 1000
        # Dimensionless frequency at which to switch to geometric optics
        geometric_optics_frequency = 1001
