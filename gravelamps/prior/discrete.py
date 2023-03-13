"""
Discrete Priors

Implements priors for discrete distributions as needed for the Morse phase which may take one of
three distinct values, or the millilensing number of images which must obviously be discrete

Written by Mick Wright
           Ania Liu
           Justin Janquart
"""

import numpy as np

from bilby.core.prior import Categorical

class UniformMorse(Categorical):
    """
    Discrete Uniform Prior for handling the Morse phase. 

    This is a restricted subset of the Categorical prior to the cases of 0, 0.5, 1

    Attributes
    ----------
    name : str
        Name of the parameter used, defaults to 'morse_phase'
    latex_label : str
        The latex compatible output to be used on plots, etc. Defaults to '$n$' 
    unit : str
        Unit of the parameter

    Methods
    -------
    rescale
        Maps the continuous distribution 0 to 1 to discrete distribution of 0, 0.5, 1
    """

    def __init__(self, name="morse_phase", latex_label="$n$", unit=None):
        super().__init__(name, 3, latex_label)
        self.categories = [0, 0.5, 1]

    def rescale(self, val):
        """
        Rescale a sample from the unit line element to one of the categories
        
        Parameters
        ----------
        val : Union[float, int, array_like]
            Uniform proabability between 0 and 1

        Returns
        -------
        Union[float, array_lik]
        """

        return np.floor(self.maximum * val)/2.
