#! /home/marcus/anaconda3/envs/netsci_cuda_11_7/bin/python
#
# Quick python script to calculate the average structure of 
# a protein (SERCA) after running an MD simulation 
# Used to generate the average structures for the pathfinding 
# analysis. 
#

import mdtraj as md
import numpy as np 

import argparse

# Most of the parameters are entered in the command line arguments 
# i.e. Stride, trajectory, topology, restart files, etc. 
def read_command_line():
    parser = argparse.ArgumentParser(description="Read in topology (parm7 or PDB) and trajectory (DCD) file. ")
    parser.add_argument('-top', type=str, required=True, help='Path to the PDB file')
    parser.add_argument('-x', type=str, required=True, help='Path to the trajectory file')
    parser.add_argument('-s', type=int, default = 1, required=False, help='Stride to read file')
    parser.add_argument('-out', type=str, required=True, help='Path to save the output PBD file')

    

    args = parser.parse_args()

    top_file = args.top
    dcd_file = args.x
    output = args.out

    if top_file[-3:] != "pdb" and top_file[-5:] != "parm7":
        print('Input file for the topology must be a .pdb or .parm7 file, not a .{} file'.format(top_file[-6:].split(".")[1]))
        print('Exiting code...')
        exit()
    if dcd_file[-3:] != "dcd" and dcd_file[-3:] != ".nc":
        print('Input file for the trajectory must be a .dcd file, not a .{} file'.format(dcd_file[-3:]))
        print('Exiting code...')
        exit()
    if output[-3:] != "pdb":
        print('Input file for the trajectory must be a .pdb file, not a .{} file'.format(output[-3:]))
        print('Appending .pdb to filepath. ')
        output += ".pdb"
        print('Continuing')
    
    return top_file, dcd_file, args.s, output

top_file, dcd_file, stride, output = read_command_line()

traj = md.load(dcd_file, top = top_file, stride = stride)

# average_coords = np.mean(traj.xyz, axis 
avg_struct = np.array([np.mean(traj.xyz, axis = 0)])

average_traj = md.load_frame(dcd_file, top = top_file, index = 0)

average_traj.xyz = avg_struct

average_traj.save(output)


