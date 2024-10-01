### This file is to run the generalized correlation
### calculation using the MI  algorithm in NETSCI.
###
### Created by Marcus Hock 09/22/2023
### Edited online 10/26/2023

import mdtraj as md # On TSCC this needs to be imported first. 
import argparse
import cuarray
import netchem
import netcalc
import numpy as np
from time import time

# Most of the parameters are entered in the command line arguments 
# i.e. Stride, trajectory, topology, restart files, etc. 
def read_command_line():
    parser = argparse.ArgumentParser(description="Read in topology (PDB) and trajectory (DCD) file. ")
    parser.add_argument('-pdb', type=str, required=True, help='Path to the PDB file')
    parser.add_argument('-dcd', type=str, required=True, help='Path to the DCD file')
    parser.add_argument('-start', type=int, required=True, help='Starting frame number (normally 0)')
    parser.add_argument('-stop', type=int, required=True, help='Last frame number (normally lenght of traj -1)')
    parser.add_argument('-stride', type=int, default=1, help='Stride to analyze the trajectory at')
    parser.add_argument('-rest_r', type=str, default=None, help='If running a restart, the path to the R restart npy')
    parser.add_argument('-rest_ab', type=str, default=None, help='If running a restart, the path to the R restart npy')
    parser.add_argument('-method', type=str, default='loc', help="Type of alignment: 'loc' or 'none'. Default is 'loc'")


    args = parser.parse_args()

    pdb_file = args.pdb
    dcd_file = args.dcd

    if pdb_file[-3:] != "pdb":
        print('Input file for the topology must be a .pdb file, not a .{} file'.format(pdb_file[-3:]))
        print('Exiting code...')
        exit()
    if dcd_file[-3:] != "dcd":
        print('Input file for the trajectory must be a .dcd file, not a .{} file'.format(dcd_file[-3:]))
        print('Exiting code...')
        exit()
    return pdb_file, dcd_file, args.start, args.stop, args.stride, args.rest_r, args.rest_ab, args.method

def localize_coordinates(trajectory_file, topology_file, protein_graph, num_frames):

    local_dist_cutoff = 1.0 # nm

    # use the coordinates from the 
    num_nodes = protein_graph.numNodes()

    locally_aligned_nodes = np.zeros((num_nodes, 3*num_frames)).astype(np.float32)
    print(f"loading files {trajectory_file} and {topology_file}.")

    # Load in first frame 
    first_frame_traj = md.load_frame(trajectory_file, top = topology_file, index = 0)

    # Select CA atoms only 
    CA_atoms = first_frame_traj.top.select('name CA')

    # Select non-protein additional atoms (only one per residue though)
    other_res_list = []
    additional_atoms = []
    for index in first_frame_traj.top.select('not protein'):
        if first_frame_traj.top.atom(index).residue not in other_res_list:
            additional_atoms.append(index)
            other_res_list.append(first_frame_traj.top.atom(index).residue)
            

    input_atoms = np.hstack((CA_atoms, additional_atoms))
    CA_only = False

    if CA_only:
        traj = md.load(trajectory_file, top=topology_file, atom_indices = input_atoms)
        print("Used the CA atoms only and the first atoms of the nucleotides")

    
    else:
        # Reload the first frame to create a new topology with CA ony 
        traj = md.load_frame(trajectory_file, index = 0, top=topology_file, atom_indices = input_atoms)

        # Using the default read in coordinates from netsci to apply to the mdtraj instance 
        COM_coordintes = protein_graph.nodeCoordinates().toNumpy2D()
        new_coords = np.zeros((num_frames, num_nodes, 3))

        temp = traj.xyz[0]
        other_temp = COM_coordintes[0,0]

        # x coords
        new_coords[:,:,0] = COM_coordintes[:,0:num_frames].T
        # y coords
        new_coords[:,:,1] = COM_coordintes[:,num_frames:2*num_frames].T
        # z coords
        new_coords[:,:,2] = COM_coordintes[:,2*num_frames:3*num_frames].T

        # Divide by ten because netsci uses Angstroms as distances, and MD traj nm

        traj.xyz=new_coords.astype(np.float32)/10

    print("constructing local alignments.")
    for i, atom1 in enumerate(traj.topology.atoms):
        first_frame_coords = traj.xyz[0,:,:]
        atom1_coords = first_frame_coords[i, :]
        close_atom_indices = []
        for j, atom2 in enumerate(traj.topology.atoms):
            if i == j: continue
            atom2_coords = first_frame_coords[j, :]
            dist = np.linalg.norm(atom2_coords - atom1_coords)
            if dist <= local_dist_cutoff:
                close_atom_indices.append(j)
    
        traj.superpose(traj, atom_indices=close_atom_indices, ref_atom_indices=close_atom_indices)
        
        
        positions = traj.xyz[:,i,:] - traj.xyz[0,i,:]
        locally_aligned_nodes[i, 0:num_frames] = positions[:,0]
        locally_aligned_nodes[i, num_frames:2*num_frames] = positions[:,1]
        locally_aligned_nodes[i, 2*num_frames:3*num_frames] = positions[:,2]    
    
    protein_graph.nodeCoordinates().fromNumpy2D((locally_aligned_nodes*10).astype(np.float32))

    return protein_graph

def main():
    start_time = time()
    # Read in file names
    topology_file, dcd_file, start_frame, stop_frame, stride, r_file, ab_file, method = read_command_line()
    
    # Need to comment out after testing 
    '''topology_file = '/home/marcus/Documents/SERCA/serca-atp-datp/simulation_structures/stripped_trajectories/replicates_from_3W5A/3W5A_original_ATP_protein_align.pdb'
    dcd_file= '/home/marcus/Documents/SERCA/serca-atp-datp/simulation_structures/stripped_trajectories/replicates_from_3W5A/3W5A_original_ATP_protein_align_50_25ns.dcd'
    start_frame = 0
    stop_frame= 100
    stride= 5
    r_file= None
    ab_file= None
    method = 'loc'''

    save_file_str = topology_file.split('/')[-1][0:-4]+"_{}_{}".format(start_frame,stop_frame)
    # Initialize graph for analysis
    print("Creating netchem graph")
    protein_graph = netchem.Network()


    protein_graph.init(
        trajectoryFile=dcd_file,
        topologyFile=topology_file,
        firstFrame=start_frame,
        lastFrame=stop_frame,
        stride=stride
    )

    # Define cuArrays for calculation
    # Correlation Matrix R
    protein_R = cuarray.FloatCuArray()
    # Define pairs to compute 
    protein_ab = cuarray.IntCuArray()

    # Count the number of nodes in the analysis (generally number of residues)
    protein_num_nodes = protein_graph.numNodes()
    protein_num_node_pairs = protein_num_nodes ** 2
    protein_ab.init(
        protein_num_node_pairs,
        2,
    )
    for i in range(protein_num_nodes):
        for j in range(protein_num_nodes):
            protein_node_pair_index = i * protein_num_nodes + j
            protein_ab[protein_node_pair_index][0] = i
            protein_ab[protein_node_pair_index][1] = j

    # Define constants for MI analysis
    # These values likely should not change
    
    # Number of data points
    protein_n = protein_graph.numFrames()
    # Dimmensionality of data
    protein_d = 3
    protein_xd = 2
    # K-nearest neighbors 
    # Can use 4. 
    # k= 6 recommended by lane 
    protein_k = 6 
    protein_platform = netcalc.GPU_PLATFORM
    if protein_platform == 0:
        print("Using CUdA (value = {})".format(protein_platform))

    # Check to see if the local alignment should be used (by default yes)
    if method == 'loc':
        protein_graph = localize_coordinates(dcd_file, topology_file, protein_graph, protein_n)

    # Check to see if there was a restart file supplied for the calcualtion
    if r_file == None or ab_file == None:
        print("No restart files used. Starting new calculation.")
        netcalc.generalizedCorrelationWithCheckpointing(
            X=protein_graph.nodeCoordinates(),
            R=protein_R,
            ab=protein_ab,
            k=protein_k,
            n=protein_n,
            d=protein_d,
            xd=protein_xd,
            platform=protein_platform,
            checkpointFrequency=1000,
            checkpointFileName=str(save_file_str + "_check"),
        )

    else:
        print("Restart AB and R files provided, running MI calculation from previous data...")
        print("Using the following R file: \n{}".format(r_file))
        print("Using the following ab file: \n{}".format(ab_file))
        netcalc.generalizedCorrelationRestartWithCheckpointing(
            X=protein_graph.nodeCoordinates(),
            R=protein_R,
            k=protein_k,
            n=protein_n,
            d=protein_d,
            xd=protein_xd,
            platform=protein_platform,
            checkpointFrequency=1000, # Picked a default value of 1000
            checkpointFileName=str(save_file_str + "_check"),
            restartAbFileName=ab_file,
            restartRFileName=r_file,
        )





    # Reshape the output from the MI calculation to a numpy NxN array for N nodes (residues)
    protein_R_np = protein_R.toNumpy2D().reshape(
        protein_num_nodes,
        protein_num_nodes,
    )
    saved_MI_name = '{}_{}_{}.npy'.format(topology_file.split('/')[-1][0:-4], start_frame, stop_frame)
    print('Saving correlation matrix to :')
    print(saved_MI_name)
    np.save(saved_MI_name, protein_R_np)
    print("Saved successfully.")
    compute_time = time() - start_time
    print("Finished in {:.3f} hours. ({} seconds)".format(compute_time/3600, compute_time))





if __name__ == '__main__':

    main()
