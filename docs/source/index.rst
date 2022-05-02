Fieldkit documentation
========================

Fieldkit is a Python package for working with field files used by OpenFTS (and PolyFTS).
The goal of fieldkit is to make analyzing, plotting, and manipulating field files easy 
and painless.

Quickstart
=============

1. Install dependencies: 

    * numpy, scipy, matplotlib (for core functionality)
    * pytest (for unit tests)
    * sphinx (for documentation)
    * mdtraj, gsd, scikit-image, pandas, numba (for some specialized functionality)
    * Paraview (for VTK visualization)

2. Update your `PYTHONPATH` to contain **fieldkit**::
  
    export PYTHONPATH=<path>/fieldkit:$PYTHONPATH

3. Check that everything works ::

    python3 -c "import fieldkit as fk; field = fk.Field(npw_Nd=(32,32))"

Take a look at the :ref:`/tutorials/00_fieldkit_basics.ipynb` tutorial to learn more.


.. toctree::
   :caption: Tutorials

   tutorials/00_fieldkit_basics
   tutorials/01_fieldkit_and_openfts
   tutorials/02_manipulate_fields
   usage_from_shell_scripts

.. toctree::
   :caption: Reference
   :maxdepth: 1

  
.. toctree::
   :maxdepth: 2
   :caption: Fieldkit API

   fieldkit

.. toctree::
   :maxdepth: 1
   :caption: Additional examples

   add_ellipse
   change_resolution
   rep_field
   roll
   create_field
   initialize_phases
   visualize
 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
