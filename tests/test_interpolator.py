'''
Gravelamps Interpolator Tests

The unittest class here tests that the interpolators generate consistent data with what came
before. They do not test this deeply, larger tests should be used when changes are made, but
they will flag immediate oddities.

Written by Mick Wright 2022
'''

import os
import subprocess
import unittest

import numpy as np

class PointTests(unittest.TestCase):
    '''
    Point Lens Interpolator Testing Class
    '''

    @classmethod
    def setUpClass(cls):
        cls.executable_name = "gravelamps_generate_lens"
        cls.point_ini = "test_configs/point_interpolator.ini"
        cls.interpolator_data_folder = "interpolator_pregenerated/point"

        cls.dimensionless_frequency_file = "w.dat"
        cls.source_position_file = "y.dat"

        cls.outdir = "pointtest"
        cls.data_folder = "pointtest/data"

        subprocess.run([cls.executable_name, cls.point_ini], check=True)

    @classmethod
    def tearDownClass(cls):
        subprocess.run(["rm", "-rf", cls.outdir], check=True)

    def test_file_existence(self):
        '''
        Checks that the files were generated as expected
        '''

        assert os.path.isdir(self.outdir)
        assert os.path.isdir(self.data_folder)

        assert os.path.isfile(f"{self.data_folder}/analysis_dimensionless_frequency.dat")
        assert os.path.isfile(f"{self.data_folder}/analysis_source_position.dat")
        assert os.path.isfile(f"{self.data_folder}/analysis_amplification_factor_real.dat")
        assert os.path.isfile(f"{self.data_folder}/analysis_amplification_factor_imag.dat")

    def test_grid_file_consistency(self):
        '''
        Checks that the grid files are what they should be
        '''

        test_dimensionless_frequency_array = np.loadtxt(self.dimensionless_frequency_file)
        test_source_position_array = np.loadtxt(self.source_position_file)

        candidate_dimensionless_frequency_array =\
            np.loadtxt(f"{self.data_folder}/analysis_dimensionless_frequency.dat")
        candidate_source_position_array =\
            np.loadtxt(f"{self.data_folder}/analysis_source_position.dat")

        assert np.array_equal(test_dimensionless_frequency_array,
                              candidate_dimensionless_frequency_array)
        assert np.array_equal(test_source_position_array,
                              candidate_source_position_array)

    def test_data_file_consistency(self):
        '''
        Checks that the data files are what they should be
        '''

        test_real_file =\
            np.loadtxt(f"{self.interpolator_data_folder}/analysis_amplification_factor_real.dat")
        test_imag_file =\
            np.loadtxt(f"{self.interpolator_data_folder}/analysis_amplification_factor_imag.dat")

        candidate_real_file =\
            np.loadtxt(f"{self.data_folder}/analysis_amplification_factor_real.dat")
        candidate_imag_file =\
            np.loadtxt(f"{self.data_folder}/analysis_amplification_factor_imag.dat")

        assert test_real_file.shape == candidate_real_file.shape
        assert test_imag_file.shape == candidate_imag_file.shape
        assert np.array_equal(test_real_file, candidate_real_file)
        assert np.array_equal(test_imag_file, candidate_imag_file)

class SISTests(PointTests):
    '''
    SIS Testing Class
    '''

    @classmethod
    def setUpClass(cls):
        cls.executable_name = "gravelamps_generate_lens"
        cls.sis_ini = "test_configs/sis_interpolator.ini"
        cls.interpolator_data_folder = "interpolator_pregenerated/sis"

        cls.dimensionless_frequency_file = "w.dat"
        cls.source_position_file = "y.dat"

        cls.outdir = "sistest"
        cls.data_folder = "sistest/data"

        subprocess.run([cls.executable_name, cls.sis_ini], check=True)

class NFWTests(PointTests):
    '''
    NFW Testing Class
    '''

    @classmethod
    def setUpClass(cls):
        cls.executable_name = "gravelamps_generate_lens"
        cls.nfw_ini = "test_configs/nfw_interpolator.ini"
        cls.interpolator_data_folder = "interpolator_pregenerated/nfw"

        cls.dimensionless_frequency_file = "w.dat"
        cls.source_position_file = "y.dat"

        cls.outdir = "nfwtest"
        cls.data_folder = "nfwtest/data"

        subprocess.run([cls.executable_name, cls.nfw_ini], check=True)

if __name__ == "__main__":
    unittest.main()
