"""
Gravelamps Installation Procedures
----------------------------------

Installs Gravelamps including C++ extensions.

Gravelamps is largely a standard python package augmented with C++ libraries that form the basis of
the lensing interpolator generation. The code for this is located in gravelamps/model with the
built libraries being in gravelamps/model/lib
"""

import setuptools

with open("README.md", "r", encoding="utf-8") as readme_contents:
    long_description = readme_contents.read()

with open("requirements.txt", "r", encoding="utf-8") as requirement_file:
    requirements = requirement_file.read().split("\n")

setuptools.setup(
    name = "Gravelamps",
    url = "https://git.ligo.org/mick.wright/Gravelamps",
    use_scm_version = True,
    setup_requires = ["setuptools_scm"],
    author = "Mick Wright, Martin Hendry",
    maintainer = "Mick Wright",
    author_email = "mick.wright@ligo.org",
    license = "MIT",
    description = ("Software package designed for running template based analysis of lensed"
                   "gravitational wave signals to determine the lens profile model. Built"
                   "on top of the parameter estimation framework, Bilby, and C arbitrary"
                   "precision library arb"),
    long_description = long_description,
    packages = [
        "gravelamps",
        "gravelamps.core"
        "gravelamps.lensing"
    ],
    install_requires = requirements,
    entry_points = {
        "console_scripts" : [
            "gravelamps_inference=gravelamps.inference:main",
            "gravelamps_generate_lens=gravelamps.generate_lens:main",
            "gravelamps_generate_interpolator_data=gravelamps.lensing.generic:main"],
        "asimov_pipelines": [
            "gravelamps=gravelamps.asimov:Gravelamps"]
        },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        ],
    python_requires = ">=3.9"
)
