import setuptools

setuptools.setup(
    name="gwlensing",
    version="0.1.0",
    author="Mick Wright",
    author_email="m.wright.3@research.gla.ac.uk",
    description="Suite of tools to simulate lensed gravitational wave signals and perform parameter estimation on them using bilby",
    long_description="LONG DESCRIPTION HERE",
    url="https://github.com/mick-wright/gwlensing",
    packages=[
        "gwlensing",
        "gwlensing.inference",
        "gwlensing.lensing"
    ],
    entry_points = {
        'console_scripts': ['gwlensing_inference=gwlensing.inference.inference:main',
        'gwlensing_inference_pipe = gwlensing.inference.inference_pipe:main', 
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: Linux",
    ],
    python_requires='>=3.8',
)
