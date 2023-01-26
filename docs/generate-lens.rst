=============
Generate Lens
=============

The Generate Lens program (``gravelamps_generate_lens``) is one of the programs available within Gravelamps. When run it will generate lens data for the scenario given to it within the INI. 

        **Scope of this document**

        This document covers the usage and running of the Generate Lens program. For INI documentation, see :doc:`ini`.

Usage of the Program
====================

The program may be run as follows::

        gravelamps_generate_lens -o --option INI

The only positional argument to the Generate Lens program is the ``INI`` argument. This should be the filepath of the INI configuration file to be used by the program.

The optional arguments are as follows:

-h, --help                      Shows the help message and exits
-i, --injection                 Will use injection lensing settings instead of analysis lensing settings
-l, --local                     Will perform the lens generatioon locally instead of only generating condor submit
-s, --submit                    Will directly submit the condor job
-v, --verbose                   Displays the maximum level of information from the logger, overriding INI settings
