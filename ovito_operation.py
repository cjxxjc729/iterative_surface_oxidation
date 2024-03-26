#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

from ovito import data
from ovito.io import import_file, export_file
from ovito.data import SurfaceMesh
from ovito.modifiers import ConstructSurfaceModifier, CreateBondsModifier, ExpandSelectionModifier,InvertSelectionModifier,DeleteSelectedModifier, ClusterAnalysisModifier,SelectTypeModifier, ExpressionSelectionModifier
import sys
import os
import shutil

from ase.io import read,write
import numpy as np


def construct_surf_clusters(file_path):

  # Load a particle set and construct the surface mesh:
  #1
  f_xyz = file_path
  pipeline0 = import_file(f_xyz)
  prefix = f_xyz.split('.xyz')[0]

  # conver to xyz with Particle Identifier
  export_file(pipeline0, f_xyz, format='xyz',columns=['Particle Identifier','Particle Type','Position.X', 'Position.Y', 'Position.Z'])  
  pipeline = import_file(f_xyz)


  #2
  pipeline.modifiers.append(ConstructSurfaceModifier(
      method = ConstructSurfaceModifier.Method.AlphaShape,
      radius = 2.9,
      identify_regions = True,
      select_surface_particles=True))

  #3
  pipeline.modifiers.append(CreateBondsModifier(
      mode=CreateBondsModifier.Mode.Pairwise))

  #set the pairwise cutoff by ase

  atoms = read(f_xyz)
  uniq_ele = np.unique(atoms.symbols)

  for i in range(len(uniq_ele)):
    for j in range(i,len(uniq_ele)):
      '''
      为了识别分子，需要利用小分子内部的成键，此处定义一下原子的成键最大值。来自ovito
      '''
      ele_i = uniq_ele[i]
      ele_j = uniq_ele[j]
      if set([ele_i, ele_j]) == {'H', 'O'}:
        #pipeline.modifiers[1].set_pairwise_cutoff('H', 'O', 1.63)
        pipeline.modifiers[1].set_pairwise_cutoff('H', 'O', 1.15) 
      elif set([ele_i, ele_j]) == {'H', 'C'}:
        pipeline.modifiers[1].set_pairwise_cutoff('H', 'C', 1.74*0.7)
      elif set([ele_i, ele_j]) == {'C', 'C'}:
        pipeline.modifiers[1].set_pairwise_cutoff('C', 'C', 2.04*0.7) 
      elif set([ele_i, ele_j]) == {'C', 'O'}:
        pipeline.modifiers[1].set_pairwise_cutoff('C', 'O', 1.932*0.7)
      elif set([ele_i, ele_j]) == {'C', 'S'}:
        pipeline.modifiers[1].set_pairwise_cutoff('C', 'S', 2.1*0.7)
      elif set([ele_i, ele_j]) == {'C', 'N'}:
        pipeline.modifiers[1].set_pairwise_cutoff('C', 'N', 1.95*0.7)
      elif set([ele_i, ele_j]) == {'H', 'H'}:
        pipeline.modifiers[1].set_pairwise_cutoff('H', 'H', 1.1*0.7)
      elif set([ele_i, ele_j]) == {'H', 'S'}:
        pipeline.modifiers[1].set_pairwise_cutoff('H', 'S', 1.8*0.7)
      elif set([ele_i, ele_j]) == {'H', 'N'}:
        pipeline.modifiers[1].set_pairwise_cutoff('H', 'N', 1.65*0.7)
      elif set([ele_i, ele_j]) == {'S', 'O'}:
        pipeline.modifiers[1].set_pairwise_cutoff('S', 'O', 1.992*0.7)
      elif set([ele_i, ele_j]) == {'S', 'N'}:
        pipeline.modifiers[1].set_pairwise_cutoff('S', 'N', 2.01*0.7)
      elif set([ele_i, ele_j]) == {'S', 'Li'}:
        pipeline.modifiers[1].set_pairwise_cutoff('S', 'Li', 2.17*0.7)
      elif set([ele_i, ele_j]) == {'O', 'O'}:
        pipeline.modifiers[1].set_pairwise_cutoff('O', 'O', 1.824*0.7)
      elif set([ele_i, ele_j]) == {'O', 'N'}:
        pipeline.modifiers[1].set_pairwise_cutoff('O', 'N', 1.842*0.7)
      elif set([ele_i, ele_j]) == {'O', 'Li'}:
        pipeline.modifiers[1].set_pairwise_cutoff('O', 'Li',2.004*0.7)
      elif set([ele_i, ele_j]) == {'N', 'N'}:
        pipeline.modifiers[1].set_pairwise_cutoff('N', 'N', 1.86*0.7)
      elif set([ele_i, ele_j]) == {'N', 'Li'}:
        pipeline.modifiers[1].set_pairwise_cutoff('N', 'Li', 2.022*0.7) 
      else:
        pipeline.modifiers[1].set_pairwise_cutoff(ele_i, ele_j, 0.7)


  pipeline.modifiers.append(ExpandSelectionModifier(
      mode = ExpandSelectionModifier.ExpansionMode.Bonded,
      iterations = 7))

  #5
  pipeline.modifiers.append(InvertSelectionModifier())
  pipeline.modifiers.append(DeleteSelectedModifier())


  #6
  pipeline.modifiers.append(ClusterAnalysisModifier(
      neighbor_mode = ClusterAnalysisModifier.NeighborMode.Bonding,
      sort_by_size = True,
      cluster_coloring = True))

  #7
  pipeline.compute()

  #8. Export the mesh to a VTK file for visualization with ParaView.
  #export_file(mesh, 'surface_mesh.vtk', 'vtk/trimesh')

  export_file(pipeline, prefix+'_adsorbates_highlight.xyz', format='xyz',columns=['Particle Identifier','Particle Type','Position.X', 'Position.Y', 'Position.Z','Cluster','Color.R','Color.G','Color.B'])

  export_file(pipeline, prefix+'_adsorbates.xyz', format='xyz',columns=['Particle Identifier','Cluster'])

  #9 rename and edit prefix+'_adsorbates.txt'
  os.rename(prefix+'_adsorbates.xyz', prefix+'_adsorbates.txt')

  f_adbtxt = prefix+'_adsorbates.txt'
  with open(f_adbtxt, 'r') as f:
    lines = f.readlines()
  lines = lines[2:]
  with open(f_adbtxt, 'w') as f:
    f.writelines(lines)


def del_surf_clusters(file_path):

  # Load a particle set and construct the surface mesh:
  #1
  f_xyz = file_path
  pipeline0 = import_file(f_xyz)
  prefix = f_xyz.split('.xyz')[0]

  # conver to xyz with Particle Identifier
  export_file(pipeline0, f_xyz, format='xyz',columns=['Particle Identifier','Particle Type','Position.X', 'Position.Y', 'Position.Z'])
  pipeline = import_file(f_xyz)


  #2
  pipeline.modifiers.append(ConstructSurfaceModifier(
      method = ConstructSurfaceModifier.Method.AlphaShape,
      radius = 2.9,
      identify_regions = True,
      select_surface_particles=True))

  #3
  pipeline.modifiers.append(CreateBondsModifier(
      mode=CreateBondsModifier.Mode.Pairwise))

  #set the pairwise cutoff by ase

  atoms = read(f_xyz)
  uniq_ele = np.unique(atoms.symbols)

  for i in range(len(uniq_ele)):
    for j in range(i,len(uniq_ele)):

      ele_i = uniq_ele[i]
      ele_j = uniq_ele[j]
      if set([ele_i, ele_j]) == {'H', 'O'}:

        pipeline.modifiers[1].set_pairwise_cutoff('H', 'O', 1.1)

      elif set([ele_i, ele_j]) == {'H', 'C'}:

        pipeline.modifiers[1].set_pairwise_cutoff('H', 'C', 1.1)

      elif set([ele_i, ele_j]) == {'C', 'C'}:

        pipeline.modifiers[1].set_pairwise_cutoff('C', 'C', 1.3)

      elif set([ele_i, ele_j]) == {'C', 'O'}:

        pipeline.modifiers[1].set_pairwise_cutoff('C', 'O', 1.4)

      else:

        pipeline.modifiers[1].set_pairwise_cutoff(ele_i, ele_j, 0.7)

  pipeline.modifiers.append(ExpandSelectionModifier(
      mode = ExpandSelectionModifier.ExpansionMode.Bonded,
      iterations = 7))

  #4
  pipeline.modifiers.append(DeleteSelectedModifier())


  #5
  pipeline.compute()


  export_file(pipeline, prefix+'_with_surfatom_delled.xyz', format='xyz',columns=['Particle Identifier','Particle Type','Position.X', 'Position.Y', 'Position.Z'])


def clear_gas_mol(file_path):

  # Load a particle set and construct the surface mesh:
  #1
  f_xyz = file_path
  pipeline = import_file(f_xyz)
  prefix = f_xyz.split('.xyz')[0]


  #2
  pipeline.modifiers.append(CreateBondsModifier(
      mode = CreateBondsModifier.Mode.Uniform,
      cutoff = 3))


  #3
  pipeline.modifiers.append(ClusterAnalysisModifier(
      neighbor_mode = ClusterAnalysisModifier.NeighborMode.Bonding,
      sort_by_size = True,
      cluster_coloring = True))


  #4
  pipeline.modifiers.append(ExpressionSelectionModifier(
      expression = 'Cluster==1'))

  #compute
  pipeline.compute()
  export_file(pipeline, prefix+'_unbonded_mol_identify.xyz', format='xyz',columns=['Particle Identifier','Particle Type','Position.X', 'Position.Y', 'Position.Z','Cluster','Color.R','Color.G','Color.B'])


  #5

  pipeline.modifiers.append(InvertSelectionModifier())
  pipeline.modifiers.append(DeleteSelectedModifier())


  #6
  pipeline.compute()

  #shutil.copy(f_xyz, prefix+'before_deletion.xyz')
  export_file(pipeline, prefix+'.xyz', format='xyz',columns=['Particle Type','Position.X', 'Position.Y', 'Position.Z'])




if __name__ == "__main__":

  file_path = sys.argv[1]
  construct_surf_clusters(file_path)
  

