import setuptools

setuptools.setup(
    name="gwlensing",
    version="0.0.2",
    author="Mick Wright",
    author_email="m.wright.3@research.gla.ac.uk",
    description="Suite of tools to set up lensing runs using bilby for various lens models",
    long_description="LONG DESCRIPTION HERE",
    url="URL_HERE",
    packages=[
        "gwlensing",
        "gwlensing.inference",
        "gwlensing.lensing"
    ],
    entry_points = {
        'console_scripts': ['gwlensing_test=gwlensing.inference.testing:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
)
