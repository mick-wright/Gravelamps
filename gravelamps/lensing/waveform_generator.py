"""Gravelamps Waveform Generators

Following contains class definitions for Gravelamps waveform generators. These are child classes of
the standard bilby WaveformGenerator class.

Written by Mick Wright
           Isaac C. F. Wong 2022

Classes
-------
LensedWaveformGenerator
    Standard total lensed waveform generator used throughout Gravelamps
"""

import importlib

import numpy as np

from bilby.core.utils import infer_parameters_from_function
from bilby.gw.waveform_generator import WaveformGenerator

from gravelamps.core.conversion import (frequency_to_dimensionless_frequency,
                                        lens_mass_to_redshifted_lens_mass)
from gravelamps.core.gravelog import gravelogger

class LensedWaveformGenerator(WaveformGenerator):
    """
    Standard total lensed waveform generator used throughout Gravelamps

    This waveform generator will seek a lensing module from the waveform arguments from which to
    extract the amplification factor calculation function it will then use during calculation of
    the strain data.

    Attributes
    ----------
    lens_module : ModuleType
        The module used to extract lensing amplification factor function from
    lens_parameters : dict
        Additional lens model parameters and their values

    See Also
    --------
    bilby.gw.waveform_generator.WaveformGenerator
        Parent class from which most methods are inhereted
    """

    def __init__(self, duration=None, sampling_frequency=None, start_time=0,
                 frequency_domain_source_model=None, time_domain_source_model=None,
                 parameters=None, parameter_conversion=None, waveform_arguments=None):

        #Perform main initialisation from parent class
        super().__init__(duration, sampling_frequency, start_time, frequency_domain_source_model,
                         time_domain_source_model, parameters, parameter_conversion,
                         waveform_arguments)

        self.lens_module = importlib.import_module(waveform_arguments["lens_module"])

        if waveform_arguments["lens_module"] == "gravelamps.lensing.interpolator":
            self.lens_module.generate_interpolator(waveform_arguments["dimensionless_frequency"],
                                                   waveform_arguments["source_position_file"],
                                                   waveform_arguments["amplification_factor_real"],
                                                   waveform_arguments["amplification_factor_imag"])

        if hasattr(self.lens_module, "load_table"):
            load_table_func = getattr(self.lens_module, "load_table")
            load_table_func(waveform_arguments["lookup_table_location"])

        if hasattr(self.lens_module, "set_scaling"):
            scaling_setter_func = getattr(self.lens_module, "set_scaling")
            scaling_setter_func(waveform_arguments["scaling_constant"])

        if hasattr(self.lens_module, "amplification_factor"):
            self.amplification_factor_func = getattr(self.lens_module, "amplification_factor")
        else:
            self.amplification_factor_func = None
            gravelogger.warning("No Amplification Factor Function detected, \
                                 signal will be unlensed")

        self.source_parameter_keys.update(self.lens_parameters())

    def _strain_from_model(self, model_data_points, model):
        unlensed_waveform = model(model_data_points, **self.parameters)

        if self.amplification_factor_func is None:
            return unlensed_waveform

        if "lens_mass" in self.source_parameter_keys:
            redshifted_lens_mass =\
                lens_mass_to_redshifted_lens_mass(self.parameters["lens_mass"],
                                                  self.parameters["lens_fractional_distance"],
                                                  self.parameters["luminosity_distance"])
            dimensionless_frequency_array =\
                frequency_to_dimensionless_frequency(model_data_points, redshifted_lens_mass)
            amplification_factor =\
                self.amplification_factor_func(dimensionless_frequency_array,
                                               self.parameters["source_position"])
        elif "k" in self.source_parameter_keys:
            image_times, luminosity_distances, phases =\
                self.lens_module.gather_parameter_lists(self.lens_parameters(), self.parameters)

            amplification_factor =\
                self.amplification_factor_func(model_data_points,
                                               int(self.parameters["k"]),
                                               image_times,
                                               luminosity_distances,
                                               phases)

        else:
            amplification_factor = self.amplification_factor_func(model_data_points,
                                                                  **self.parameters)

        lensed_waveform = {}
        for key, value in unlensed_waveform.items():
            lensed_waveform[key] = np.multiply(value, amplification_factor)
        return lensed_waveform

    @property
    def lens_parameters(self):
        """
        Additional lens model parameters and their value

        Returns
        -------
        dict
            Contains the parameters as keys and the associated values

        See Also
        --------
        Each model should contain instructions on the parameters required.
        """

        return self._lens_parameters

    def _lens_parameters(self):
        if hasattr(self.lens_module, "_lens_parameters"):
            lens_parameters = self.lens_module._lens_parameters
        elif hasattr(self.lens_module, "get_lens_parameters"):
            parameter_func = getattr(self.lens_module, "get_lens_parameters")
            lens_parameters = parameter_func(self.waveform_arguments)
        else:
            lens_parameters = infer_parameters_from_function(self.amplification_factor_func)
        return lens_parameters
