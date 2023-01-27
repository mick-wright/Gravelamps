====================
Gravelamps Inference
====================

The Gravelamps Inference program (``gravelamps_inference``) is the main program within Gravelamps. When run, it will generate all the necessary lensing data using the other programs available within Gravelamps before starting a nested sampling parameter estimation run using the ``bilby`` framework. Multiple runs with the same data but differing lens models may be performed in order to perform model selection between the differing lens models.

        **Scope of this document**

        This document covers the usage and running of the Inference program. For INI documentation, see :doc:`ini`.

Usage of the Program
====================

The program may be run as follows::

        gravelamps_inference -o --option INI

The only positional argument to the Inference program is the ``INI`` argument. This should be the filepath of the INI configuration file to be used by the program. 

The optional arguments are as follows:

-h, --help               Shows the help message and exits
-i, --injection          Run contains an injection, will generate injection data
-l, --local              Will perform the run locally instead of only generating the condor submit files
-s, --submit             Will submit the condor jobs directly instead of relying on the user to do this
-v, --verbose            Displays the maximum level of information from the logger, overriding INI settings
