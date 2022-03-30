Getting Started with fieldkit
=============================

Dependencies
------------

The following Python packages and software must be installed:

   * numpy - to run most fieldkit functions
   * matplotlib or gnuplot - to model 1D fields
   * Paraview - to model 2D and 3D fields (prefer to install via conda?)
   * Pytest (optional) - to run unit tests in the tests folder

Installing Paraview
-------------------

**On Ubuntu/Debian**:: 

   sudo apt-get install paraview

Setting PYTHONPATH
------------------

If Python does not recognize fieldkit as an installed Python package

* For one session::

   export PYTHONPATH="$HOME/FTS-tools"

* By default for any terminal session::

   vim ~/.bashrc
   export PYTHONPATH="$HOME/FTS-tool"

