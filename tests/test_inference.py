'''
Gravelamps Inference Testing

Here is the unittest class and functions designed to test the inference program

Written by Mick Wright 2022
'''

import os
import subprocess
import unittest

class InferenceTests(unittest.TestCase):
    '''
    Inference Testing Class
    '''

    @classmethod
    def setUpClass(cls):
        cls.executable_name = "gravelamps_inference"
        cls.injection_ini = "test_configs/injection_test.ini"
        cls.event_ini = "test_configs/event_test.ini"

        os.mkdir("event_data")
        with open("event_data/H1.dat", "w", encoding="utf-8"):
            pass
        with open("event_data/L1.dat", "w", encoding="utf-8"):
            pass
        with open("event_data/V1.dat", "w", encoding="utf-8"):
            pass

    @classmethod
    def tearDownClass(cls):
        subprocess.run(["rm", "-rf", "injection_test"], check=True)
        subprocess.run(["rm", "-rf", "event_test"], check=True)
        subprocess.run(["rm", "-rf", "event_data"], check=True)

    def test_event(self):
        '''
        Runs a test run using an event example config
        '''
        subprocess.run([self.executable_name, self.event_ini], check=True)
        assert os.path.isfile("event_test/submit/gravelamps_inference.dag")

    def test_injection(self):
        '''
        Runs a test run using an event injection config
        '''
        subprocess.run([self.executable_name, self.injection_ini], check=True)
        assert os.path.isfile("injection_test/submit/gravelamps_inference.dag")

if __name__ == "__main__":
    unittest.main()
