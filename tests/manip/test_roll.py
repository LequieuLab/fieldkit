import fieldkit as fk

#1d
#npw = [64]
#field = fk.Field(npw_Nd=npw)
#fk.add_ellipse(field,center=[0.5], axis_lengths=[0.25])
#fk.write_to_file("fields.in",[field])

#2d
#npw = (64,64)
#field = fk.Field(npw_Nd=npw)
#fk.add_ellipse(field,center=(0.0,0.0), axis_lengths=(0.3,0.2),height=1)
#fk.write_to_file("fields.in",[field])

# 2d multiple ellipse
npw = (64,64)
field = fk.Field(npw_Nd=npw)
fk.add_ellipse(field,center=(0.5,0.5), axis_lengths=(0.3,0.2),height=1)
fk.write_to_file("fields_init.dat",[field])
fields_new = fk.roll([field],shift=(0.25,0.5))
fk.write_to_file("fields_roll.dat",fields_new)
