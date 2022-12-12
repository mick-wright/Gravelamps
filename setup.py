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
        model_subfolder = f"{current_directory}/gravelamps/model"

        #Change to the subdirectory, run the makefile, return to present
        os.chdir(model_subfolder)
        subprocess.run(['make'], check=True)
        os.chdir(current_directory)

        #Run the remaining build process
        build_ext.run(self)

with open("README.md", "r", encoding="utf-8") as readme_contents:
    long_description = readme_contents.read()

with open("requirements.txt", "r", encoding="utf-8") as requirement_file:
    requirements = requirement_file.read().split("\n")

setuptools.setup(
    name = "Gravelamps",
    url = "https://arxiv.org/abs/2112.07012",
    download_url = "https://git.ligo.org/michael.wright/Gravelamps",
    use_scm_version = True,
    setup_requires = ['setuptools_scm'],
    author = "Mick Wright, Martin Hendry",
    maintainer = "Mick Wright",
    author_email = "mick.wright@ligo.org",
    license = "MIT",
    description = ("Software package designed for running template based analysis of lensed"
                   "gravitational wave signals to determine the lens profile model. Built"
                   "on top of the parameter estimation framework, Bilby, and arbitrary"
                   "precision library arb"),
    long_description = long_description,
    packages = [
        "gravelamps",
        "gravelamps.core",
        "gravelamps.lensing"
    ],
    has_ext_modules = lambda: True,
    install_requires = requirements,
    cmdclass = {
        "build_ext": Build,
    },
    entry_points = {
        "console_scripts": [
            "gravelamps_inference=gravelamps.inference:main",
            "gravelamps_generate_lens=gravelamps.generate_lens:main",
            "gravelamps_generate_interpolator_data=gravelamps.lensing.generic:main"],
        "asimov.pipelines": ["gravelamps = gravelamps.asimov.asimov:Gravelamps"]},
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Operating System :: POSIX :: Linux",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: GW Lensing Community",
        "Topic :: GW Lensing :: Lens Mass Profile Model Selection"
        ],
    python_requires = ">=3.8",
)
