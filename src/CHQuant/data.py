import os, glob
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from .utils import Trajectory

def df_from_path(inputpath): 
    """
    This function takes the path to an XYZ trajectory or a directory of XYZs 
    and outputs a dataframe of convex hull data
    (number of points, area/volume sums/averages, atomic areas/volumes, percentage contributions to total area/volume)
    Inputs:
        inputpath (string): path to an XYZ file or directory of XYZ files
    Outputs:
        df (pandas dataframe) containing convex hull data
    """
    if os.path.isdir(inputpath):
        xyz_paths = glob.glob(f"{inputpath}/*.xyz")
    if os.path.isfile(inputpath):
        xyz_paths = [inputpath]
    
    column_volumes = []
    volume_sums = []
    volume_avgs = []
    point_counts = []
    atomic_areas = []
    area_avgs = []
    area_sums = []
    avg_densities = []
    
  
    atomic_volume_contributions_column = []
    atomic_area_contributions_column = []

    for i,path in enumerate(xyz_paths):
        
        #Reading xyz
        traj_coords, atom_names = Trajectory(path)
        
        #Convex hull, atomic properties
        point_counts.append(int(len(traj_coords)))
        volumes = []
        areas = []
        for atomindex in range(len(atom_names)): #loop over atoms 
            selected_coords = [i[atomindex] for i in traj_coords] #selecting one atom
            hull = ConvexHull(selected_coords)
    #         print(f"{atom_names[atomindex]}: {hull.volume}")
            volumes.append(hull.volume)
#             print(hull.volume)
            areas.append(hull.area)
        atomic_areas.append(areas)
        area_sums.append(sum(areas))
        area_avgs.append(sum(areas)/len(areas))
        volume_sums.append(sum(volumes))
        volume_avgs.append(sum(volumes)/len(volumes))
        column_volumes.append(volumes)
        avg_densities.append(len(volumes)/sum(volumes))
        
        #Atomic % contribution to total volume, area
        atomic_volume_contributions = []
        atomic_area_contributions = []
        
        for atomindex in range(len(atom_names)): #loop over atoms 

            selected_coords = [c[atomindex] for c in traj_coords] #selecting one atom
            hull = ConvexHull(selected_coords)
    #         print(f"{atom_names[atomindex]}: {hull.volume}")

            #Calculate percentage contributions from each atom
            atomic_volume_contributions.append(100*hull.volume/volume_sums[i])
            atomic_area_contributions.append(100*hull.area/area_sums[i])

        atomic_volume_contributions_column.append(atomic_volume_contributions)
        atomic_area_contributions_column.append(atomic_area_contributions)
        
    #data for each atom, i.e. Volumes and surface areas

    data = [volumes+atomic_areas[i]+atomic_volume_contributions_column[i]+atomic_area_contributions_column[i] for i,volumes in enumerate(column_volumes)]

    column_labels = [os.path.basename(i)[:-4] for i in xyz_paths]
    
    #Ordering of index labels should match ordering of atomic data
    index_labels = [i+" volume" for i in [atom_names][0]]+[i+" area" for i in [atom_names][0]]+[i+" % volume" for i in [atom_names][0]]+[i+" % area" for i in [atom_names][0]]

    other_data = [avg_densities,volume_avgs, volume_sums,area_avgs,area_sums,point_counts]
    other_data_labels = ["Density","Volume avg","Volume sum","Area avg", "Area sum","Number of points"]

    for i in other_data_labels:
        index_labels.insert(0,i)

    for i in range(len(data)): #Loop over number of trajectories loaded
        for stat in other_data:
            data[i].insert(0,stat[i])

    data = np.array(data).T

    df = pd.DataFrame(data, index=index_labels,columns=column_labels)

    return df  
    
def df_from_points(points,system_name="system_name",atom_names=None): 
    """
    This function takes an array of points and outputs a dataframe of convex hull data of that point cloud
    (number of points, area/volume sums/averages, atomic areas/volumes, percentage contributions to total area/volume);
    unlike df_from_path which reads the system and atom names from the input filename,
    you also need to specify the system name and atom names for this function
    Inputs:
        points (array): NxFx3 array of coordinates where N is the number of points and F is the number of frames in the trajectory
        ^this should now match the output of load_directory_convex
        system_name (string): name of the system; will be used as the column title of the dataframe
        atom_names (list): list of atom names i.e. C1, O2... to use for the row labels
    Outputs:
        df (pandas dataframe) containing convex hull data
    """
    column_volumes = []
    volume_sums = []
    volume_avgs = []
    point_counts = []
    atomic_areas = []
    area_avgs = []
    area_sums = []
    avg_densities = []
    i=0 #since this is adapted from the old df function
  
    atomic_volume_contributions_column = []
    atomic_area_contributions_column = []

    #Convex hull, atomic properties
    
    volumes = []
    areas = []

    for atomindex in range(len(atom_names)): #loop over atoms 
    
        #Needs fixing
        #Way to select coordinates depends on input array shape? could make separate functions; comment out as appropriate
        # selected_coords = [i[atomindex] for i in points] #selecting one atom 
        selected_coords = points[atomindex]
        point_counts.append(int(len(selected_coords)))
        # selected_coords = points[atomindex]

        #print(selected_coords)
        # print(len(selected_coords))
        hull = ConvexHull(selected_coords)
#         print(f"{atom_names[atomindex]}: {hull.volume}")
        volumes.append(hull.volume)
#             print(hull.volume)
        areas.append(hull.area)
    atomic_areas.append(areas)
    area_sums.append(sum(areas))
    area_avgs.append(sum(areas)/len(areas))
    volume_sums.append(sum(volumes))
    volume_avgs.append(sum(volumes)/len(volumes))
    column_volumes.append(volumes)
    avg_densities.append(len(volumes)/sum(volumes))

    #Atomic % contribution to total volume, area
    atomic_volume_contributions = []
    atomic_area_contributions = []

    for atomindex in range(len(atom_names)): #loop over atoms 

        selected_coords = points[atomindex]
        hull = ConvexHull(selected_coords)
#         print(f"{atom_names[atomindex]}: {hull.volume}")

        #Calculate percentage contributions from each atom
        atomic_volume_contributions.append(100*hull.volume/volume_sums[i])
        atomic_area_contributions.append(100*hull.area/area_sums[i])

    atomic_volume_contributions_column.append(atomic_volume_contributions)
    atomic_area_contributions_column.append(atomic_area_contributions)
        
    #data for each atom, i.e. Volumes and surface areas

    data = [volumes+atomic_areas[i]+atomic_volume_contributions_column[i]+atomic_area_contributions_column[i] for i,volumes in enumerate(column_volumes)]

#     column_labels = [os.path.basename(i)[:-4] for i in xyz_paths]
    column_labels = [system_name]
    
    #Ordering of index labels should match ordering of atomic data
    index_labels = [i+" volume" for i in [atom_names][0]]+[i+" area" for i in [atom_names][0]]+[i+" % volume" for i in [atom_names][0]]+[i+" % area" for i in [atom_names][0]]

    other_data = [avg_densities,volume_avgs, volume_sums,area_avgs,area_sums,point_counts]
    other_data_labels = ["Density","Volume avg","Volume sum","Area avg", "Area sum","Number of atoms"]

    for i in other_data_labels:
        index_labels.insert(0,i)

    for i in range(len(data)): #Loop over number of trajectories loaded
        for stat in other_data:
            data[i].insert(0,stat[i])

    data = np.array(data).T

    df = pd.DataFrame(data, index=index_labels,columns=column_labels)

    return df

    
def df_from_coords(points,system_name="system_name",atom_names=None): 
    """
    This function takes an array of points and outputs a dataframe of convex hull data of that trajectory.
    (number of points, area/volume sums/averages, atomic areas/volumes, percentage contributions to total area/volume);
    unlike df_from_path which reads the system and atom names from the input filename,
    you also need to specify the system name and atom names for this function
    Inputs:
        points (array): NxFx3 array of coordinates where N is the number of points and F is the number of frames in the trajectory
        ^this should now match the output of load_directory_convex
        system_name (string): name of the system; will be used as the column title of the dataframe
        atom_names (list): list of atom names i.e. C1, O2... to use for the row labels
    Outputs:
        df (pandas dataframe) containing convex hull data
    """
    column_volumes = []
    volume_sums = []
    volume_avgs = []
    point_counts = []
    atomic_areas = []
    area_avgs = []
    area_sums = []
    avg_densities = []
    i=0 #since this is adapted from the old df function
  
    atomic_volume_contributions_column = []
    atomic_area_contributions_column = []

    #Convex hull, atomic properties
    
    volumes = []
    areas = []

    for atomindex in range(len(atom_names)): #loop over atoms 
    
        #Needs fixing
        #Way to select coordinates depends on input array shape? could make separate functions; comment out as appropriate
        selected_coords = [i[atomindex] for i in points] #selecting one atom 
        # selected_coords = points[atomindex]
        point_counts.append(int(len(selected_coords)))
        # selected_coords = points[atomindex]

        #print(selected_coords)
        # print(len(selected_coords))
        hull = ConvexHull(selected_coords)
#         print(f"{atom_names[atomindex]}: {hull.volume}")
        volumes.append(hull.volume)
#             print(hull.volume)
        areas.append(hull.area)
    atomic_areas.append(areas)
    area_sums.append(sum(areas))
    area_avgs.append(sum(areas)/len(areas))
    volume_sums.append(sum(volumes))
    volume_avgs.append(sum(volumes)/len(volumes))
    column_volumes.append(volumes)
    avg_densities.append(len(volumes)/sum(volumes))

    #Atomic % contribution to total volume, area
    atomic_volume_contributions = []
    atomic_area_contributions = []

    for atomindex in range(len(atom_names)): #loop over atoms 

        # selected_coords = points[atomindex]
        selected_coords = [i[atomindex] for i in points]
        hull = ConvexHull(selected_coords)
#         print(f"{atom_names[atomindex]}: {hull.volume}")

        #Calculate percentage contributions from each atom
        atomic_volume_contributions.append(100*hull.volume/volume_sums[i])
        atomic_area_contributions.append(100*hull.area/area_sums[i])

    atomic_volume_contributions_column.append(atomic_volume_contributions)
    atomic_area_contributions_column.append(atomic_area_contributions)
        
    #data for each atom, i.e. Volumes and surface areas

    data = [volumes+atomic_areas[i]+atomic_volume_contributions_column[i]+atomic_area_contributions_column[i] for i,volumes in enumerate(column_volumes)]

#     column_labels = [os.path.basename(i)[:-4] for i in xyz_paths]
    column_labels = [system_name]
    
    #Ordering of index labels should match ordering of atomic data
    index_labels = [i+" volume" for i in [atom_names][0]]+[i+" area" for i in [atom_names][0]]+[i+" % volume" for i in [atom_names][0]]+[i+" % area" for i in [atom_names][0]]

    other_data = [avg_densities,volume_avgs, volume_sums,area_avgs,area_sums,point_counts]
    other_data_labels = ["Density","Volume avg","Volume sum","Area avg", "Area sum","Number of points"]

    for i in other_data_labels:
        index_labels.insert(0,i)

    for i in range(len(data)): #Loop over number of trajectories loaded
        for stat in other_data:
            data[i].insert(0,stat[i])

    data = np.array(data).T

    df = pd.DataFrame(data, index=index_labels,columns=column_labels)

    return df