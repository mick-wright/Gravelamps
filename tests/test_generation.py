'''
Gravelamps Lens Generation Testing

Here is the unittest class and functions designed to test the lens generation program

Written by Mick Wright 2022
'''

import os
import subprocess
import unittest

import numpy as np

class LensGenerationTests(unittest.TestCase):
    '''
    Lens Generation Testing Class
    '''

    @classmethod
    def setUpClass(cls):
        cls.executable_name = "gravelamps_generate_lens"
        cls.cluster_generation_ini = "test_configs/cluster_generation.ini"
        cls.local_generation_ini = "test_configs/local_generation.ini"

    @classmethod
    def tearDownClass(cls):
        subprocess.run(["rm", "-rf", "cluster_generation"], check=True)
        subprocess.run(["rm", "-rf", "local_generation"], check=True)

    def test_cluster_generation(self):
        '''
        Performs a test run using cluster generation
        '''
        subprocess.run([self.executable_name, self.cluster_generation_ini], check=True)
        assert os.path.isfile("cluster_generation/submit/generate_analysis_interpolator_data.sub")

    def test_local_generation(self):
        '''
        Performs a test run using local generation
        '''
        subprocess.run([self.executable_name, self.local_generation_ini], check=True)

        dimensionless_frequency_file = "local_generation/data/analysis_dimensionless_frequency.dat"
        source_position_file = "local_generation/data/analysis_source_position.dat"
        amplification_factor_real_file =\
            "local_generation/data/analysis_amplification_factor_real.dat"
        amplification_factor_imag_file =\
            "local_generation/data/analysis_amplification_factor_imag.dat"

        assert os.path.isfile(dimensionless_frequency_file)
        assert os.path.isfile(source_position_file)
        assert os.path.isfile(amplification_factor_real_file)
        assert os.path.isfile(amplification_factor_imag_file)

        dimensionless_frequency_array = np.loadtxt(dimensionless_frequency_file)
        source_position_array = np.loadtxt(source_position_file)
        amplification_factor_real_array = np.loadtxt(amplification_factor_real_file)
        amplification_factor_imag_array = np.loadtxt(amplification_factor_imag_file)

        assert amplification_factor_real_array.shape == (len(source_position_array),
                                                         len(dimensionless_frequency_array))
        assert amplification_factor_imag_array.shape == amplification_factor_real_array.shape

if __name__ == "__main__":
    unittest.main()
