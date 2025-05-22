#!/usr/bin/env python3

import argparse
import fieldkit as fk
import numpy as np

parser = argparse.ArgumentParser(description='calculate the dense and dilute phase densities from density data')
parser.add_argument('-i', '--input', type=str, help='filename for density data')
parser.add_argument('-s', '--step_size', type=int, default=3, help='stepsize for calculating the first derivative, if data is noisey try increasing this')
parser.add_argument('--dense_trim', type=int, default=0, help='amount to trim the dense phase')
parser.add_argument('--dilute_trim', type=int, default=0, help='amount to trim the dilute phase')
parser.add_argument('--plot', action='store_true', help='plot the results')
parser.add_argument('--field_index', type=int, default=0, help='field index to use')

args = parser.parse_args()
input_filename = args.input
step_size = args.step_size
dense_trim = args.dense_trim
dilute_trim = args.dilute_trim
plot = args.plot
field_index = args.field_index

try:  # Handle the output from an OpenFTS simulation.
    density_fields = fk.read_from_file(input_filename)
    centered_density_fields = fk.center_1D_density_fields(density_fields)
    density_field = centered_density_fields[field_index]
except:  # Handle the output from an MD simulation processed by particletools.
    density_trajectory = np.load(input_filename)
    density_field_trajectory = []
    nframes = density_trajectory.shape[0]
    for frame in range(nframes // 2, nframes):  # Average over the last 50%.
        density_data = density_trajectory[frame]
        density_field = fk.Field(npw=[density_data.shape[0]],
                                 data=density_data[:, 1],
                                 h=np.asarray([[density_data[-1, 0] -\
                                                density_data[0, 0]]]))
        density_field.coords = density_data[:, 0]
        density_field_trajectory.append(density_field)
    
    # For noisey data, like the profiles of frames, using a larger step_size
    # helps when centering each frame.

    centered_density_fields = fk.center_1D_density_fields(density_field_trajectory,
                                                          step_size=10)

    # Create a density field which is the average of the centered frames.

    density_data = 0
    for centered_density_field in centered_density_fields:
        density_data += centered_density_field.data.real
    density_data /= len(centered_density_fields)
    density_field = centered_density_fields[0]
    density_field.data.real = density_data

# Find the phase boundaries and write the phase data.

output_filename = f"{input_filename[:input_filename.find('.')]}.phasedat"
fk.find_phase_densities(density_field, step_size, dense_trim, dilute_trim, 
                        plot, output_filename)

