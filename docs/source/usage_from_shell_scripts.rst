Usage from shell scripts
==========================

Fieldkit also contains a number of stand-alone scripts that can be called from shell scripts. 
These can be useful for to performing common operations without needing to separately boot up Python.

All of these scripts are located in the `<path-to-fieldkit>/tools/` directory. If you plan on using these scripts frequently, it is recommended that you update your `PATH`::

    export PATH=<path-to-fieldkit>/tools/

All of these scripts have names beginning with `fk_*`. Here's one example converting a `density.dat` file to a VTK file. ::

    <path>/tools/fk_convert_to_vtk.py density.dat

Take a look around the `tools/` directory to see the other functionalities that have been implemented.
