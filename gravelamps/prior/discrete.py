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

class ImageNumber(Categorical):
    """
    Discrete Uniform Prior for handling the number of millilensing image signals.

    Lightly modified such that the minimum must be 1

    Attributes
    ---------
    ncategories : int
        Number of potential images
    name : str
        Name of the parameter used, defaults to 'num_images'
    latex_label : str
        The latex compatible output to be used on plots, etc. Defaults to '$n_{\mathrm{signals}}$'.
    unit : str
        Unit of the parameter
    """ 

    def __init__(self, ncategories, name="num_images", 
                 latex_label="$n_{\mathrm{signals}}$", unit=None):
        super().__init__(ncategories, name=name, latex_label=latex_label, unit=unit)
        maximum = ncategories + 1
        minimum = 1
        self.categories = np.arange(self.minimum, self.maximum)


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
        super().__init__(3, name=name, latex_label=latex_label)
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
