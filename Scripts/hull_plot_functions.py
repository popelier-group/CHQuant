import os,glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


from scipy.spatial import ConvexHull
import plotly.graph_objects as go
from sklearn.preprocessing import StandardScaler

#Ichor functions
from ichor.core.files.xyz import XYZ, Trajectory
from ichor.core.calculators.connectivity.distance_connectivity_calculator import *
from ichor.core.calculators.features.alf_features_calculator import calculate_alf_features
from ichor.core.calculators.alf import get_atom_alf, calculate_alf_atom_sequence
from ichor.core.calculators import calculate_alf_cahn_ingold_prelog

from ichor.core.calculators.c_matrix_calculator import calculate_c_matrix
from ichor.core.atoms.alf import ALF 
from ichor.core.atoms.atoms import Atoms, Atom
# from ichor.core.calculators import *

#Bienfait concave hull script
# from concave_hull_script import estimate_volume, estimate_surface_area

def ms(x, y, z, radius, resolution): #generates sphere vertices for the atoms
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)
    
def plot_loaded_convex(fig,atom_names, point_clouds, hull_vertices): #plot hulls after loading from directory
    #shows all CHs
        
    colours = {"C":"black", "N":"blue", "H":"dimgray","O":"red"}
    
    selected_atom_types = ["C"]
    

    
    for atom, points, vertices in zip(atom_names, point_clouds, hull_vertices):
        atomtype = atom[0] #assuming one-letter element name
        # print(atomtype)
        fig.add_trace(go.Mesh3d(x=vertices[:, 0], 
                            y=vertices[:, 1], 
                            z=vertices[:, 2], 
                            color=colours[atomtype], 
                            opacity=.2, #originally 0.5, then 0.2
                            lighting=dict(facenormalsepsilon=0,ambient=0.5, specular=1.0,fresnel=0.5),
                            alphahull=0))


            #Mist plot of all points; change allpoints to vertices to plot vertices

        fig.add_trace(go.Scatter3d(
        x=points[:, 0], 
        y=points[:, 1], 
        z=points[:, 2],mode='markers',showlegend=False,
 #             name=f"{traj.atom_names[atomindex]}",
            marker=dict(color=colours[atomtype],size=0.5)))
        fig.update_layout(  #hide axes
        scene = dict(
            xaxis = dict(visible=False),
            yaxis = dict(visible=False),
            zaxis = dict(visible=False),),)
            
    return fig
    
def plot_mist(fig,atom_names, point_clouds, hull_vertices): #plot hulls after loading from directory
    #shows all CHs
        
    colours = {"C":"black", "N":"blue", "H":"dimgray","O":"red"}
    
    selected_atom_types = ["C"]
    

    
    for atom, points, vertices in zip(atom_names, point_clouds, hull_vertices):
        atomtype = atom[0] #assuming one-letter element name
        # print(atomtype)
        


            #Mist plot of all points; change allpoints to vertices to plot vertices

        fig.add_trace(go.Scatter3d(
        x=points[:, 0], 
        y=points[:, 1], 
        z=points[:, 2],mode='markers',showlegend=False,
 #             name=f"{traj.atom_names[atomindex]}",
            marker=dict(color=colours[atomtype],size=0.5)))
        fig.update_layout(  #hide axes
        scene = dict(
            xaxis = dict(visible=False),
            yaxis = dict(visible=False),
            zaxis = dict(visible=False),),)
            
    return fig
    
def add_molecule_spheres(fig,path,frame,radius,resolution, bondwidth=20):
    #Molecule using spheres instead of markers
    #Draw molecule
    traj = Trajectory(path)
    traj_coords = traj.coordinates
    atomlist = [i for i in traj.atom_names]

    colours = {"C":"black", "N":"blue", "H":"dimgray","O":"red"}

    geometry = traj_coords[frame]
    atom_names = traj.atom_names

    fig = fig

    for coords,atom_name in zip(geometry,atom_names):
        atomtype=atom_name[0]
        x,y,z = ms(coords[0],coords[1],coords[2],radius,resolution)
        fig.add_trace(go.Surface(x=x,y=y,z=z, opacity=1,colorbar=None, #colorbar is property (need to load in)
        colorscale=[[0, colours[atomtype]], [1, colours[atomtype]]],
        name=f"{traj.atom_names[int(atom_name[1:])-1]}",
        showlegend=True,
        #Comment out this line for no lighting
        lighting=dict(specular=0.5,ambient=0.5,fresnel=0.5)
        ))
        
#         fig.update_traces(showscale=False)
        fig.update_traces()
    #Could use these to get circles in the legend
#     for coords,atom_name in zip(geometry,atom_names):
#         atomtype = atom_name[0]
#         fig.add_trace(go.Scatter3d(x=np.array(coords[0]),y=np.array(coords[1]),z=np.array(coords[2]),visible="legendonly",
#         mode='markers',name=f"{traj.atom_names[int(atom_name[1:])-1]}",marker=dict(color=colours[atomtype],size=15)))

    connectivity = connectivity_calculator_distance(traj[0])
    connectivity = np.triu(connectivity)

    columns, rows = np.where(connectivity==1) #pairs of connected atoms
    pairs = [i for i in zip(columns,rows)]

    for pair in pairs:
        atom1 = geometry[pair[0]]
        atom2 = geometry[pair[1]]
        coords = np.array([atom1,atom2])
        fig.add_trace(go.Scatter3d(x=coords[:, 0],y=coords[:, 1],z=coords[:, 2],mode="lines",marker=dict(size=0,color=None),showlegend=False,line=dict(color="gray",width=bondwidth)))
    
#     fig.update_layout(
#         scene = dict(
#             xaxis = dict( range=[-8,8],visible=False),
#             yaxis = dict( range=[-8,8],visible=False),
#             zaxis = dict( range=[-8,8],visible=False),),)
    
    #Format legend
#     fig.update_layout(legend=dict(itemwidth=50, # change it to 0.3
#     #                               entrywidthmode='fraction',
#                                   orientation='h',
#     #                               y=1.2,
#                                   xanchor="right", yanchor="middle",
#                                   x=0.4))
    return fig