import os, glob
import numpy as np
from scipy.spatial import ConvexHull

from .utils import Trajectory

def generate_atomic_convex_hulls(xyzpath): 
    """
    This function takes an XYZ trajectory and generates convex hulls for each atom in the molecule
    Inputs:
        xyzpath (string): path to an XYZ trajectory
    Outputs:
        hulls (list): a list of scipy ConvexHull objects
    """
    traj_coords, atom_names = Trajectory(xyzpath)
  
    hulls = []
    for atomindex,atom_name in enumerate(atom_names):
        selected_coords = [i[atomindex] for i in traj_coords] #get coordinates of just one atom in the trajectory
        hull = ConvexHull(selected_coords)
        hulls.append(hull)
    return hulls
    
def save_vertices_points(hull,filename,outputdir): 
    """
    This function takes a scipy ConvexHull object and saves its point cloud and convex hull vertices' coordinates as CSVs
    Inputs:
        hull (scipy ConvexHull object)
        filename (string): name that you want to save the file as, e.g. an atom name like C1, O2, H3 etc.
        outputdir (string): directory to save the CSVs into
    Outputs:
        filename_points.csv (CSV): coordinates of the point cloud
        filename_convex_vertices.csv (CSV): coordinates of the convex hull vertices 
    """
    vertices = hull.points[hull.vertices]
    np.savetxt(f"{outputdir}/{filename}_convex_vertices.csv", vertices, delimiter=",") 
    np.savetxt(f"{outputdir}/{filename}_points.csv", hull.points, delimiter=",")
    
def load_directory_convex(path): 
    """
    This function takes a directory with the outputs from save_vertices_points and loads the point clouds and convex hull vertices
    Inputs:
        path (string): directory containing the CSV outputs from save_vertices_points
    Outputs:
        atom_names (list): the filenames of the CSVs that were loaded; these would be the atom names if the hulls were saved using load_traj.save_points()
        point_clouds (list): list of arrays of point clouds; each array is Nx3 where N is the number of points
        hull_vertices (list): list of arrays of vertices; each array is Nx3 where N is the number of vertices
    """
    hull_vertices = []
    point_clouds = []
    
    points_paths = glob.glob(f"{path}/*points.csv")
    vertices_paths = glob.glob(f"{path}/*convex_vertices.csv")
    
    points_paths = sorted(points_paths,key=lambda x:int(os.path.basename(x).split("_")[0][1:])) #sort based on atom index
    vertices_paths = sorted(vertices_paths,key=lambda x:int(os.path.basename(x).split("_")[0][1:]))
    
    for i,j in zip(points_paths,vertices_paths): #Check that ordering is the same
        assert os.path.basename(i).split("_")[0]==os.path.basename(j).split("_")[0]
        
    atom_names = [os.path.basename(i).split("_")[0] for i in points_paths]
#     print(atom_names)
    
    for atom, points, vertices in zip(atom_names,points_paths,vertices_paths):
        point_clouds.append(np.genfromtxt(points, delimiter=','))        
        hull_vertices.append(np.genfromtxt(vertices, delimiter=','))  
        
    return atom_names, point_clouds, hull_vertices
    #point_clouds should have 20k points each; vertices should have fewer points
      