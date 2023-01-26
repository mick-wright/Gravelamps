"""
Gravelamps Core
---------------

The Gravelamps core module handles the internal backend functions used to create and run the main
program functionalitis within Gravelamps. These are not expected to be necessary on the user side
to be necessary for the running of individual analyses. They are exposed to the user for the
purpose of allowing modification of core fucntionality to extend the utility of Gravelamps and the
development of additional models in a uniform fashion

Mick Wright 2022

Submodules
----------
conversion
    Handles conversion of the units of physical quantities
file_handling
    Handles input and output files to Gravelamps programs
gravelog
    Controls the implemtntation of logging within Gravelamps programs
graveparser
    Controls the parsing of commandline arguments to Gravelamps programs
module_handling
    Handles the loading and checking of modules in an agnostic sense
"""

from . import conversion
from . import file_handling
from . import gravelog
from . import graveparser
from . import module_handling
