Initializing 3D Phases with 2 or 3 species
==========================================

.. automodule:: fieldkit.field_initialize
   :members: initialize_phase, add_ellipse
   :undoc-members:
   :show-inheritance:

Examples
--------

The following example outputs vtk files for Field objects created by the initialize_phase function.::

   import numpy as np

   #parameters
   npw = (32,32,32)
   h = np.eye(3) * 2

   new_field = fk.initialize_phase("A15", npw, h)
   fk.write_to_VTK("A15.vtk", new_field)

   new_field = fk.initialize_phase("alt-C15", npw, h)
   fk.write_to_VTK("alt-C15.vtk", new_field)a

This examples shows how to add one of more ellipses in 1D/2D.::

   #1d
   npw = [64]
   field = fk.Field(npw_Nd=npw)
   fk.add_ellipse(field,center=[0.5], axis_lengths=[0.25])
   fk.write_to_file("fields.in",[field])

   #2d
   npw = (64,64)
   field = fk.Field(npw_Nd=npw)
   fk.add_ellipse(field,center=(0.0,0.0), axis_lengths=(0.3,0.2),height=1)
   fk.write_to_file("fields.in",[field])

   # 2d multiple ellipse
   npw = (64,64)
   field = fk.Field(npw_Nd=npw)
   fk.add_ellipse(field,center=(0.5,0.5), axis_lengths=(0.3,0.2),height=1)

   fk.add_ellipse(field,center=(0.4,0.5), axis_lengths=(0.05,0.05),height=-1,smoothing_width=0.02)
   fk.add_ellipse(field,center=(0.6,0.5), axis_lengths=(0.05,0.05),height=-1,smoothing_width=0.02)
   fk.write_to_file("fields.in",[field])
