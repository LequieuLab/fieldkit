# fieldkit
fieldkit is a Python package used to visualize and manipulate 1D, 2D, and 3D fields. 

## Table of Contents
1. [Getting Started](#get_started)
	1. [Prerequisites](#prereqs)
	2. [Installing Paraview](#install_paraview)
	3. [Setting PYTHONPATH](#set_python)
2. [Creating a Field through fieldkit](#create_field)
	1. [Field class](#field_class)
	2. [Reading and Writing Functions](#read_write_functions)
		1. [read_from_file](#read_from_file)
		2. [write_to_file](#write_to_file)
		3. [write_to_VTK](#write_to_VTK)
		4. [Example 1](#example1)
3. [Initializing 3D Phases with 2 or 3 species](#intialize)
	1. [intialize_phase](#intialize_phase)
	2. [Example 2](#example2)
4. [Manipulating Fields](#manipulate)
	1. [change_resolution](#change_resolution)
	2. [replicate_fields](#replicate_fields)
	3. [Example 3](#example3)
	4. [Example 4](#example4)
5. [Visualizing Fields](#visualize)
	1. [1D Field Using matplotlib or gnuplot](#plot)
	2. [2D/3D Fields Using paraview](#paraview)

## Getting Started <a name="get_started"></a>
### Prequisites <a name="prereqs"></a>
Python 3 and above is required to use this package

The following Python packages and software must be installed:
* numpy
* matplotlib or gnuplot - model 1D fields
* Paraview - model 2D and 3D fields

### Installing Paraview <a name="install_paraview"></a>
On Ubuntu/Debian:
`sudo apt-get install paraview`

### Setting PYTHONPATH <a name="set_python"></a>
If Python does not recognize fieldkit as an installed Python package

###### For one session 
`export PYTHONPATH="$HOME/FTS-tools"`

###### By default for any terminal session
```
vim ~/.bashrc
export PYTHONPATH="$HOME/FTS-tool"
```

## Creating a Field through fieldkit <a name="create_field"></a>
### Field class <a name="field_class"></a>
Class for Field object used to model 1D, 2D, 3D structures.  

##### Attributes:
* npw_Nd: number of grid points in each dimension
* npw_total: total number of grid points for all dimensions
* dim: dimension of field
* data: stored values of field at each grid point
* h: a cell tensor, in which the columns are box vectors that describe the coordinate system
* is_real_space: tells if Fourier transform has been used
* coords: stores x,y,z coordinates of each grid point
   
### Reading and Writing Functions<a name="read_write_functions"></a>
All functions are found  in field_io.py 

#### read_from_file(filename)<a name="read_from_file"></a>
> Reads from textfile and output a list of Field objects.
>
> Adapted from iotools.py.
> 
> Args:
> * filename: String for name of file to be read through this function.
>
> Returns:
> * field_list: A list of Field objects.


#### write_to_file(filename, fields)<a name="write_to_file"></a>
> Creates a text file based on the fields variable.
>	
> Adapted from iotools.py.
>  
> Args:
> * filename: String for name of file to be written through this function.
> * fields: A list of Fields objects.
>          
> Returns:
>
> 	A text file based on fields.
>  
> Raises:
> * TypeError: if fields is not a list, than make it a one element list.

#### write_to_VTK(filename, fields)<a name="write_to_VTK"></a>
> Creates a VTK file based on the fields variable.
>
> Adapted from FTS-tools/plot/PolyFTS_to_VTK.py.
>
> Args:
> * filename: String for name of file to be written through this function.
> * fields: A list of Fields objects.
>  
> Returns:
>
> 	A VTK file based on fields.
>       
> Raises: 
>
>	If VTK file already exist, the function exits 

#### Example 1<a name="example1"></a>
```
import fieldkit as fk

filename = "density.dat"

fields = fk.read_from_file(filename)
fk.write_to_file("tmp.dat", fields)
fk.write_to_VTK("density.vtk", fields)
#outputs density.vtk

```

All files are in the manual-tests/demo folder.

## Initializing 3D Phases with 2 or 3 species <a name="initialize"></a>
The following phases are implemented in the initialize_phase function.
* 2 phases
	* A15
	* C15
	* single-diamond
	* double-gyroid
* 3 phases
	* alt-C15
	* alt-diamond

This function is found in field_initialize.py

#### intialize_phase(phase, npw, h) <a name="intialize_phase"></a>
> Creates a list of Field objects for 3D phases with 2 or 3 species. 
>    
> Adapted from FTS-tools/field/init_fields_3d_2spec.py and from FTS-tools/fields/ini_fields_3d_3spec.py.
>    
> Args:
> * phase: name of the phase (see description above for phases implemented in this function)
> * npw = number of grid points for each dimension
> * h = a cell tensor, in which the columns are box vectors that describe the coordinate system
>    
> Returns:
> * field_list: list of Field objects, each object for each species.
>
> Raises:
> * NotImplementedError: if phase is not a phases implemented in this function (see description above), the function will exit.


#### Example 2 <a name="example2"></a>
Refer to <a name="example1">Example 1</a> for defination of fields.

```
import numpy as np

#parameters
npw = (32,32,32)
h = np.eye(3) * 2

new_field = fk.intialize_phase("A15", npw, h)	
fk.write_to_VTK("A15.vtk", new_field)

new_field = fk.intialize_phase("alt-C15", npw, h)	
fk.write_to_VTK("alt-C15.vtk", new_field)

#outputs A15.vtk and alt-C15.vtk
```


## Manipulating Fields <a name="manipulate"></a>
All function are found in field_manip.py

#### change_resolution(field_old, resolution_new) <a name="change_resolution"></a>
> For a list of Field objects, change the resolution of each Field object.
>
> Args:
> * field_old: a list of Field objects
> * resolution_new: a tuple that defines the new resolution of each Field object.
>   
> Returns:
> * field_new: a list of Field objects with each object set to a resolution of resolution_new.


#### replicate_field(fields, nreplicates) <a name="replicate_field"></a>
> For a list of Fields, replicate each Field object by nreplicates.
>
> Adapted from FTS-tools/replicate_fields.py and FTS-tools/lib/fieldtools.py.
>    
> Args:
> * fields: a list of Field objects
> * nreplicated: number of replicates
> Returns: 
> * fields_list: a list of Field objects, in which each Field object is replicated by nreplicates amount of times.


#### Example 3<a name="example3"></a>
Refer to <a name="example1">Example 1</a> for defination of fields.

```
field_new = fk.change_resolution(fields, (32,32))
fk.write_to_file("res_32.dat", field_new)
fk.write_to_VTK("res_32.vtk", field_new)
```

#### Example 4 <a name="example4"></a>

```
field_new = fk.replicate_fields(fields, (2,2))
fk.write_to_VTK("repfields.vtk", field_new)
```


## Visualizing Fields <a name="visualize"></a> 
### 1D Field Using matplotlib or gnuplot <a name="plot"></a>

###### matplotlib example
```
filename = "density_1D.dat"
fields = fk.ReadFromFile(filename)
newfields = fk.ReplicateFields(fields, 3)

#very basic plotting
plt.figure
plt.plot(newfields[0].coords, newfields[0].data)
plt.show()
```

###### gnuplot example
On Terminal (Linux):
```gnuplot```

```p './density_1D.dat' u 1:2```

Note: 1:2 means plotting column 1 and column 2

### 2D/3D Fields Using Paraview <a name="paraview"></a>

#### For 2D plots
Using density.vtk as an example. First, change the current directory to the same directory density.vtk is in.

On Terminal: 
```paraview density.vtk```

On Paraview:
Click "Apply" and change Representation to "Surface"


#### For 3D plots
Using A15.vtk as an example. First, change the current directory to the same directory A15.vtk is in.

On Terminal:
```paraview A15.vtk```

On Paraview:
Click "Apply" and change Representation to "Volume"


