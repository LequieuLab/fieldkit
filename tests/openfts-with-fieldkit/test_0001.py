
import openfts
import fieldkit as fk 

fts = openfts.OpenFTS()

fts.cell(cell_scale=3.8,cell_lengths=[1.0, 0.866],tilt_factors=[0.5],dimension=2,length_unit='Rg')
fts.driver(dt=5.0,nsteps=1000,output_freq=50,type='SCFT',stop_tolerance=1e-05,)
fts.field_updater(type='EMPEC',update_order='simultaneous',adaptive_timestep=False,lambdas=[1.0, 1.0])
fts.variable_cell(lambda_=1.0,stop_tolerance=0.0001,update_freq=10,shape_constraint='none')
fts.field_layout(npw=[64, 64],random_seed=1)
fts.model(Nref=1.0,bref=1.0,chiN=20.0,type='MeltChiAB')

field_filename = "test.in"
field = [fk.Field(npw_Nd = (64,64)) for i in range(2)]
fk.write_to_file(field_filename, field)
fts.init_fields_from_file(file=field_filename, field_type='species')

fts.add_molecule(alpha=1.0,ds=0.01,nblocks=2,block_fractions=[0.25, 0.75],block_species=['A', 'B'],type='PolymerContinuous',volume_frac=1.0)

fts.add_operator(averaging_type='none',type='Hamiltonian')
fts.add_operator(averaging_type='none',type='CellStress')

fts.output_default()

fts.add_species(label='A',stat_segment_length=1.0)
fts.add_species(label='B',stat_segment_length=1.0)

fts.write_to_file('input.json')

if __name__ == '__main__':
  fts.run()


# for pytest
def test_all(tmpdir):
  import os
  from pytest import approx

  # extra work to find path to fields.in
  mypath = os.getenv('PYTEST_CURRENT_TEST')
  myfile = os.path.abspath(__file__) # abs path to this test file
  mypath = '/'.join(myfile.split('/')[:-1]) # get path of test, exclude final arg (which is filename)
  import shutil
  old = mypath + '/%s' % field_filename
  new = tmpdir + '/%s' % field_filename
  shutil.copy(old,new)
  

  os.chdir(tmpdir) # go into tmp directory
  fts.run()

  H = fts.io_get_scalar_operator("H.real")
  assert(H == approx(-1.389307,abs=1e-4))

  h00 = fts.io_get_scalar_operator("h[0][0]")
  h10 = fts.io_get_scalar_operator("h[1][0]")
  h01 = fts.io_get_scalar_operator("h[0][1]")
  h11 = fts.io_get_scalar_operator("h[1][1]")
  Lx = (h00**2 + h10**2)**0.5
  Ly = (h01**2 + h11**2)**0.5

  assert(Lx == approx(3.879078,abs=1e-3))
  assert(Ly == approx(Lx,abs=1e-3))
