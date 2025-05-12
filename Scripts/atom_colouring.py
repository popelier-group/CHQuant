import numpy as np
import pandas as pd
from matplotlib import colormaps
from sklearn.preprocessing import MinMaxScaler
import plotly.graph_objects as go
from pathlib import Path
import os, glob

from ichor.core.calculators.connectivity.distance_connectivity_calculator import connectivity_calculator_distance
from ichor.core.files.xyz import XYZ, Trajectory

def ms(x, y, z, radius, resolution): #generates sphere vertices for the atoms
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)

def matplotlib_to_plotly(cmap, pl_entries): #convert matplotlib colormap to list of colors
        h = 1.0/(pl_entries-1)
        pl_colorscale = []

        for k in range(pl_entries):
            C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
            pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

        return pl_colorscale
        
'''
Required Inputs:
    fig: plotly Figure object, can be created using "fig = go.Figure()"
    path: path to the molecule's geometry/trajectory (.xyz file)
    volumes: list of property values (e.g. volume, energy, error), one value per atom
    
Options:
    colormap_name: name of colormap (https://matplotlib.org/stable/users/explain/colors/colormaps.html), default "jet"
    frame: index of the frame to use (if path is a trajectory), default 0
    radius: radius of atom spheres, default 0.2
    resolution: controls how many points are used to draw each sphere, default 20
        
Outputs:
    fig: plotly Figure object
'''

def add_molecule_spheres_volumescale(fig,path,volumes,colormap_name="jet",frame=0,radius=0.2,resolution=20,data_min=0, data_max=0,labelled=False,bondwidth=10,colorbarwidth=30,colorbarfontsize=20,colorbarfractions=False,colorbarlength=0.5):
    
    
    #Load geometry
    traj = Trajectory(path)
    traj_coords = traj.coordinates
    geometry = traj_coords[frame]
    atom_names = traj.atom_names
    assert len(volumes)==len(atom_names)

    #Scale the volumes (or other property) to be between 0 and 1
#     scaler = MinMaxScaler(feature_range=(0,1))
#     scaled_volumes = scaler.fit_transform(volumes.reshape(-1,1))
#     scaled_volumes = [i[0] for i in scaled_volumes]
    
    #Assuming the input has already been scaled since we want to use the same colorbar for all figures
    scaled_volumes = volumes
    
    #Try rounding the values if there are errors
#     scaled_volumes = [round(i,5) for i in scaled_volumes]
#     print(scaled_volumes)
    
    #Draw atoms 
    cm = colormaps[colormap_name]
    rgb_values = []
    
    annotations = []
    
    for coords,atom_name,volume in zip(geometry,atom_names,scaled_volumes):
        atomtype=atom_name[0]
        x,y,z = ms(coords[0],coords[1],coords[2],radius,resolution)
        currentcolor = f"rgb{cm(volume)[:-1]}" #Get RGB values
#         print(currentcolor)
        rgb_values.append(cm(volume))
        
        #Draw atom
        fig.add_trace(go.Surface(x=x,y=y,z=z, opacity=1,colorbar=None,
        colorscale=[[0, currentcolor], [1, currentcolor]],
        showscale=False, 
        name=f"{traj.atom_names[int(atom_name[1:])-1]}",
        showlegend=False,
        text=volume,
                                 
        #Comment this line out for no lighting
        # lighting=dict(specular=0.5,ambient=0.5,fresnel=0.5)  
        
                                )) 
        
        #Uncomment these lines to label each atom with its corresponding (original) value 
        if labelled:
            test_scaler = MinMaxScaler(feature_range=(0,1))
            test_scaler.fit_transform(np.array([data_min,data_max]).reshape(-1,1))        
            annotations.append(dict(x=coords[0], y=coords[1], z=coords[2],text=round(test_scaler.inverse_transform(np.array(volume).reshape(1,-1))[0][0],2)))
    
    #Uncomment this to get labels for each point
    fig.update_layout(scene=dict(annotations=annotations))
    
    #Uncomment this to output csv of rgb values for each atom (for doing the figures in Blender)
    rgb_df = pd.DataFrame(rgb_values, columns=["R", "G", "B","Alpha"])
    rgb_df.to_csv("rgb_values.csv")
    
    #Calculate connectivity (bonds) using ICHOR
    connectivity = connectivity_calculator_distance(traj[0])
    connectivity = np.triu(connectivity)
    columns, rows = np.where(connectivity==1) #pairs of connected atoms
    pairs = [i for i in zip(columns,rows)]

    #Draw bonds
    for pair in pairs:
        atom1 = geometry[pair[0]]
        atom2 = geometry[pair[1]]
        coords = np.array([atom1,atom2])
        fig.add_trace(go.Scatter3d(x=coords[:, 0],y=coords[:, 1],z=coords[:, 2],mode="lines",marker=dict(size=0,color=None),showlegend=False,line=dict(color="gray",width=bondwidth)))
    
    #Adding custom colorbar 
    color_range = matplotlib_to_plotly(cm, 255)
    data_range = data_max-data_min
    
    #Colorbar without zero labelled; default colorbar
    if not colorbarfractions:
        colorbar_trace  = go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         marker=dict(
                             colorscale=color_range, 
                             showscale=True,
                             cmin=-5,
                             cmax=5,
                             # colorbar=dict(thickness=10, tickvals=[-5, 5], ticktext=[round(min(volumes),2), round(max(volumes),2)], outlinewidth=0)
                             colorbar=dict(thickness=colorbarwidth, tickvals=[-5, 5], ticktext=[f"{round(data_min,2):.2f}", f"{round(data_max,2):.2f}"], outlinewidth=0, tickfont=dict(size=colorbarfontsize))
                         ),
                         hoverinfo='none'
                        )
                    
    #Colorbar with zero labelled (if you have +ve and -ve values)
    # zero_fraction = (0-data_min)/data_range #assuming the data has positive and negative values
    # zero_scaled = -5+10*zero_fraction #-5 is the min value of the colorbar, 10 is the range
    # colorbar_trace  = go.Scatter(x=[None],
                             # y=[None],
                             # mode='markers',
                             # marker=dict(
                                 # colorscale=color_range, 
                                 # showscale=True,
                                 # cmin=-5,
                                 # cmax=5,
                                 ## Original colorbar
                                 ## colorbar=dict(thickness=10, tickvals=[-5, 5],
                                               ## ticktext=[round(data_min,2), round(data_max,2)],outlinewidth=0)
                                 
                                 # Colorbar with zero labelled; need to calculate where it should be earlier                                 
                                 # colorbar=dict(thickness=10, tickvals=[-5,zero_scaled, 5],
                                 # ticktext=[round(data_min,2),"0", round(data_max,2)],outlinewidth=0)
                                 
                             # ),
                             # hoverinfo='none'
                            # )
    
    #Colorbar with intervals labelled; 9 for set1 colormap, 8 labels needed
    
    # zero_fraction = (0-data_min)/data_range #assuming the data has positive and negative values
    # zero_scaled = -5+10*zero_fraction #-5 is the min value of the colorbar, 10 is the range
    
    # colorbar_trace  = go.Scatter(x=[None],
                             # y=[None],
                             # mode='markers',
                             # marker=dict(
                                 # colorscale=color_range, 
                                 # showscale=True,
                                 # cmin=-5,
                                 # cmax=5,
                                 ## Original colorbar
                                 ## colorbar=dict(thickness=10, tickvals=[-5, 5],
                                               ## ticktext=[round(data_min,2), round(data_max,2)],outlinewidth=0)
                                 
                                 # Colorbar with zero labelled; need to calculate where it should be earlier                                 
                                 # colorbar=dict(thickness=10, tickvals=[-5,zero_scaled, 5],
                                 # ticktext=[round(data_min,2),"0", round(data_max,2)],outlinewidth=0)
                                 
                             # ),
                             # hoverinfo='none'
                            # )
    if colorbarfractions:
    
        color_count=9
        fraction = (data_max-data_min)/(color_count)
        interval_labels = [data_min+fraction*i for i in range(1,color_count+1)]
        interval_labels = [round(i,2) for i in interval_labels]
        interval_labels = [f"{i:.2f}" for i in interval_labels]
        
        fractions = [-5+i*(10/(color_count)) for i in range(color_count+1)]
        colorbar_trace  = go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         marker=dict(
                             colorscale=color_range, 
                             showscale=True,
                             cmin=-5,
                             cmax=5,                              
                             colorbar=dict(thickness=colorbarwidth, tickvals=fractions, len=colorbarlength,
                             ticktext=[round(data_min,2)]+interval_labels+ [round(data_max,2)],outlinewidth=0,tickfont=dict(size=colorbarfontsize))
                             
                         ),
                         hoverinfo='none'
                        )
    
        
    fig.add_trace(colorbar_trace)
    
    #Hide axes and legend, and set the background to white
    fig.update_layout(  
    scene = dict(
        zaxis = dict(visible=False),
        xaxis = dict(visible=False),
        yaxis = dict(visible=False)
    ),        
    showlegend=False,
    template='simple_white'
    )

    fig.update_xaxes(visible=False, showgrid=False, zeroline=False,showticklabels=False,color='rgb(0,0,0)')
    fig.update_yaxes(visible=False, showgrid=False,zeroline=False,showticklabels=False,color='rgb(0,0,0)')
    
    return fig
    
def add_molecule_spheres(fig,path,frame=0,radius=0.2,resolution=20):
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
        # print(atom_name)
        atomtype=atom_name[0] #hardcoded for single-letter atoms
        x,y,z = ms(coords[0],coords[1],coords[2],radius,resolution)
        fig.add_trace(go.Surface(x=x,y=y,z=z, opacity=1,colorbar=None, #colorbar is property (need to load in)
        colorscale=[[0, colours[atomtype]], [1, colours[atomtype]]],
        name=f"{traj.atom_names[int(atom_name[1:])-1]}",
        showlegend=False,
        lighting=dict(specular=0.5,ambient=0.5,fresnel=0.5)))
        fig.update_traces(showscale=False)
        # fig.update_traces()
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
        fig.add_trace(go.Scatter3d(x=coords[:, 0],y=coords[:, 1],z=coords[:, 2],mode="lines",marker=dict(size=0,color=None),showlegend=False,line=dict(color="gray",width=10)))
    
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
