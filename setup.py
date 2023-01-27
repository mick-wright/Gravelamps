"""Gravelamps Installation Procedures.

Following code is the installation procedure for Gravelamps. It is largely a standard setuptools
style installation, with a single exception. That being that the C++ libraries that form the basis
of the lens generation codes will be compiled during installation from a makefile procedure. This
makefile is located in gravelamps/models/Makefile.

Notes
-----

This procedure is able to fetch all python dependencies, but is unable to ascertain if the C++
dependencies are satisfied and will fail if they are not due to the failure of the Makefile. Please
make sure that these dependencies are fulfilled. A list of them can be found in the README
associated with the repository.

Written by Mick Wright 2021.

"""

import subprocess
import os

import setuptools

from setuptools.command.build_ext import build_ext

class Build(build_ext):
    """
    Build procedure.

    The standard build procedure has been modified so as to run the Makefile for the C++ libraries
    that form the basis of the lens generation code.

    Methods
    -------
    run
        Executes Makefile and runs standard build procedure

    """

    def run(self):
        """
        Executes Makefile and runs standard build procedure.

        The Makefile for the C++ libraries is located in the repository's `gravelamps/model`
        folder. This Makefile generates the C++ libraries for each of the models within. This
        procedure is unable to ascertain whether the C++ dependencies are satisfied so this must
        be done before hand otherwise installation will fail.

        """

        current_directory = os.getcwd()
        model_subfolder = f"{current_directory}/gravelamps/model"

        os.chdir(model_subfolder)
        subprocess.run(['make'], check=True)
        os.chdir(current_directory)

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
    package_data = {'gravelamps': ['model/lib/*.so']},
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
        "asimov.pipelines": ["gravelamps=gravelamps.asimov:Gravelamps"]},
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        ],
    python_requires = ">=3.8",
)
