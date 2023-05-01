============
Installation
============

Conda Instllation
=================

The simplest means of installing Gravelamps is via the conda python package manager, specifically in the conda-forge channel. To install Gravelamps in this manner from a conda enviornment run::

        $ conda install -c conda-forge gravelamps

replacing ``conda`` with ``mamba`` if using that implementation. 

Manual Installation
===================

Manual installation whilst more difficult has been written to be as easily done as possible. Firstly, clone the repository using the following command::

        $ git clone git@git.ligo.org:mick.wright/Gravelamps.git

From this point, ensure that you have all of the necessary dependencies to build Gravelamps. These dependencies are:

        #. ``GNU make``
        #. ``g++``
        #. ``arb``
        #. ``flint``
        #. ``python``
        #. ``pip``

The next step is to install all of the Gravelamps python dependencies. This can be done simply by running::

        $ pip install -r requirements.txt

from within the Gravelamps repository folder. With the requirements available, installation of Gravelamps is simply::

        $ pip install .

from within the Gravelamps repository folder. If you are installing with the intent of developing, an installation that will track changes that are made without the need to reinstall can be done by running::

        $ pip install -e .
