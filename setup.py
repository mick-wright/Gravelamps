'''
Gravelamps Installation Procedure

Mick Wright 2021
'''

import subprocess
import os

import setuptools

from setuptools.command.build_ext import build_ext

class Build(build_ext):
    '''Customised version of the build - additionally runs the makefile in the lensing subfolder'''
    def run(self):
        #Get the current directory, from which find the lensing subdirectory
        current_directory = os.getcwd()
        lensing_subfolder = current_directory + "/gravelamps/lensing"

        #Change to the subdirectory, run the makefile, return to present
        os.chdir(lensing_subfolder)
        subprocess.run(["make"], check=True)
        os.chdir(current_directory)

        #Run the remaining build process
        build_ext.run(self)

with open("README.md", "r") as readme_contents:
    long_description = readme_contents.read()

setuptools.setup(
    name = "Gravelamps",
    version = "1.0.0",
    author = "Mick Wright",
    author_email = "michael.wright@ligo.org",
    license = "MIT",
    description = ("Software package designed for running template based analysis of lensed"
                   "gravitational wave signals to determine the lens profile model. Built"
                   "on top of the parameter estimation framework, Bilby, and arbitrary"
                   "precision library arb"),
    long_description = long_description,
    url = "https://github.com/mick-wright/Gravelamps",
    packages = [
        "gravelamps",
        "gravelamps.inference",
        "gravelamps.lensing"
    ],
    has_ext_modules = lambda: True,
    cmdclass = {
        "build_ext": Build,
    },
    entry_points = {
        "console_scripts": ["gravelamps_local_inference=gravelamps.inference.inference:main",
    "gravelamps_pipe_inference=gravelamps.inference.inference_pipe:main",
    "gravelamps_generate_lens_local=gravelamps.lensing.generate_lens_local:main",
    "gravelamps_generate_lens_pipe=gravelamps.lensing.generate_lens_pipe:main"],
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Operating System :: POSIX :: Linux",
        "License :: OSI Approved :: MIT License",
        ],
    python_requires = ">=3.8",
)
