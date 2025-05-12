import os, glob
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from .data import df_from_path, df_from_points
from .utils import Trajectory
from .vertices import generate_atomic_convex_hulls, save_vertices_points, load_directory_convex
from .functional_group_list import molecule_atomcounts, molecule_functional_groups
from .functional_groups import functional_group_overlapping_volume, functional_group_nonoverlapping_volume, functional_group_overlap_check

class load_traj:
    def __init__(self, inputpath, name):
        assert os.path.isfile(inputpath) and os.path.basename(inputpath[-4:])==".xyz"
        self.name = name
        self.trajpath = inputpath
        self.traj_coords, self.atom_names = Trajectory(self.trajpath)
        self.atomic_hulls = generate_atomic_convex_hulls(self.trajpath)
        self.data_df = df_from_path(inputpath) #df_from_path can also work on directories of XYZs
        self.molecule_name = self.name.split("_")[0].lower() #assuming name is {molecule}_{method}
        self.functional_groups = molecule_functional_groups[self.molecule_name]
        
    
    def save_points(self, outputpath):
        """
        This function outputs the atom point clouds and convex hull vertices as CSVs
        Inputs:
            outputpath (string): path to a directory
        Outputs:
            atomname_points.csv (CSV): coordinates of a specific atom 
            atomname_convex_vertices.csv (CSV): coordinates of a specific atom's convex hull vertices
        """
        for hull,filename in zip(self.atomic_hulls,self.atom_names):
            save_vertices_points(hull,f"{filename}",outputpath)
    
    def functional_group_df(self):
        """
        This function outputs a dataframe containing the convex hull volumes with and without overlap for the molecule's functional groups
        Functional groups are defined in functional_group_list.py
        Inputs: 
            none
        Outputs:
            df (pandas dataframe)
        """
        rows = []
        for indices, name in zip(*molecule_functional_groups[self.name]):
            new_row = {}
            
            volume_overlap = functional_group_overlapping_volume(self.trajpath, indices)
            volume_no_overlap = functional_group_nonoverlapping_volume(self.trajpath, indices)
                        
            new_row.update({"Name":name})
            new_row.update({"Volume with overlap (volume when pooling points)":volume_overlap})
            new_row.update({"Volume with no overlap (sum of atomic volumes)":volume_no_overlap})
            rows.append(new_row)
        df = pd.DataFrame(rows)
        return df
        
class load_vertices:
    def __init__(self, inputpath, name):
        assert os.path.isdir(inputpath)
        self.name = name
        self.molecule_name = self.name.split("_")[0].lower() #assuming name is {molecule}_{method}
        self.dirpath = inputpath
        self.atom_names, self.point_clouds, self.hull_vertices = load_directory_convex(inputpath)
        
        #Very small differences when using vertices vs points, approx 10e-15, probably due to numerical instability
        self.data_df = df_from_points(self.hull_vertices, self.name, self.atom_names)
        
        
    def functional_group_atom_overlap(self):
        """
        This function outputs a dictionary with functional group names and
        lists of atom pairs with overlapping convex hulls; if the list is empty it means there are no overlaps
        Inputs:
            none
        Outputs:
            overlaps (dictionary): dictionary containing functional group names as keys and lists of tuples of pairs of atoms with overlapping convex hulls
        """
        overlaps = {}
        for indices, name in zip(*molecule_functional_groups[self.molecule_name]):            
            overlap_output = functional_group_overlap_check(self.dirpath,indices)[1]
            overlaps.update({name:overlap_output})
            # overlaps.update({name.lower():overlap_output}) #lowercase to make it easier to look up names
        return overlaps