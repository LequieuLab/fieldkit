Adding an Ellipse to an existing field
======================================

.. automodule:: fieldkit.field_initialize
:members: initialize_phase
:undoc-members:
:show-inheritance:

Example
-------

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
