import numpy as np
    
def Trajectory(xyzpath):
    """
    This function takes an XYZ file and reads the coordinates and atom names
    Inputs:
        xyzpath (string): path to an XYZ trajectory
    Outputs:
        traj_coords (numpy array): NxAx3 array where N is the number of frames, A is the number of atoms in the molecule
        atom_names (list): list of atom names i.e. C1, O2, H3...
    """
    separated_xyzs = []
    with open(xyzpath) as file:
        lines = [line.rstrip() for line in file]
        atomcount = int(lines[0])
        counter = atomcount+2 #since xyz format has 2 lines at the top
        
        #split trajectory into frames
        for y in range(len(lines)): 
            if (y+1)%counter==0:
                separated_xyzs.append(lines[y-(counter-1):y+1]) 

    traj_coords = []
    
    #split frames into coordinates
    for x in range(len(separated_xyzs)): 
        traj_coords.append(np.array([row.split()[1:] for row in separated_xyzs[x][2:]],dtype="float")) 

    traj_coords = np.array(traj_coords)
    
    atoms = [row.split()[0] for row in separated_xyzs[0][2:]] #read atom names from first frame
    atom_names = [f"{atom}{i+1}" for i,atom in enumerate(atoms)]
    
    return traj_coords, atom_names