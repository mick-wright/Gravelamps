"""
Gravelamps
----------
Gravelamps: Gravitational Wave Lensing Mass Profile Model Selection

Gravelamps is a package built upon the bilby framework designed to perform template based analyses
of both simulated and real gravitationally lensed gravitational wave signals to determine the mass
density profile of the lensing object. In so doing, it will also give estimates of both the lens
and source parameters for each of the mass density profiles that is tested.

It is able to do this in both the wave optics only regime as well as in a hybrid environment
crossing into the geometric optics regime at a specified threshold allowing a great deal of
flexibility in the mass spectrum of the lensing object. The particularly complex calculations
required to compute the wave optics case are doing using the C precision library arb with the C++
libraries that have been written using this contained within gravelamps.model.

The source code of Gravelamps as well as documentation and installation instructions may be found
at https://git.ligo.org/mick.wright/gravelamps

Lead Developer: Mick Wright (mick.wright@ligo.org)

Contributors: Martin Hendry,
              Ania Liu,
              Isaac C. F. Wong,
              Eungwang Seo
"""

from pkg_resources import get_distribution, DistributionNotFound

import asimov

from .core import conversion, file_handling, gravelog, graveparser, module_handling
from . import lensing
from . import inference

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "Development"
