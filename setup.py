import setuptools

setuptools.setup(
    name="Gravelamps",
    version="0.4.0",
    author="Mick Wright",
    author_email="m.wright.3@research.gla.ac.uk",
    description="Suite of tools to simulate lensed gravitational wave signals and perform parameter estimation on them using bilby",
    long_description="LONG DESCRIPTION HERE",
    url="https://github.com/mick-wright/gwlensing",
    packages=[
        "gravelamps",
        "gravelamps.inference",
        "gravelamps.lensing"
    ],
    entry_points={
        'console_scripts': ['gravelamps_local_inference=gravelamps.inference.inference:main',
            'gravelamps_pipe_inference=gravelamps.inference.inference_pipe:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: Linux",
    ],
    python_requires='>=3.8',
)
