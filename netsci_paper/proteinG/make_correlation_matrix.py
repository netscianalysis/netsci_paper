"""
Make the 2D generalized correlation matrices for the B1 domain
of protein G, showing correlation between residues.
"""
import time

import numpy as np
import mdtraj as md

import cuarray
import netchem
import netcalc

local_dist_cutoff = 0.75 # in nm

def residue_com(traj, res, frame=0):
    first_frame_coords = traj.xyz[frame,:,:]
    com = np.array([0.0, 0.0, 0.0])
    total_mass = 0.0
    for k, atom1 in enumerate(res.atoms):
        mass = atom1.element.mass
        com += mass * first_frame_coords[atom1.index, :]
        total_mass += mass
        
    com /= total_mass
    return com

def load_trajectory(topology_file, trajectory_file, silent=False):
    if not silent:
        print(f"loading files {trajectory_file} and {topology_file}.")
    traj = md.load(trajectory_file, top=topology_file)
    #traj = md.load(trajectory_file, top=topology_file, stride=100)
    first_frame = 0
    last_frame = traj.n_frames-1
    num_frames = last_frame - first_frame + 1
    if not silent:
        print("last_frame:", last_frame)
        print("num_frames:", num_frames)
    return traj, first_frame, last_frame, num_frames

def create_netchem_graph(topology_file, trajectory_file, first_frame, last_frame, silent=False):
    if not silent:
        print("Creating netchem graph")
    graph = netchem.Graph()
    graph.init(
        trajectoryFile=trajectory_file,
        topologyFile=topology_file,
        firstFrame=first_frame,
        lastFrame=last_frame,
    )
    num_nodes = graph.numNodes()
    num_node_pairs = num_nodes**2
    return graph, num_nodes, num_node_pairs

def locally_align_trajectory(traj, first_frame, last_frame, num_frames, num_nodes, silent=False):
    locally_aligned_nodes = np.zeros((num_nodes, 3*num_frames)).astype(np.float32)
    if not silent:
        print("constructing local alignments.")
    for i, res1 in enumerate(traj.topology.residues):
        atom1_coords = residue_com(traj, res1)
        close_atom_indices = []
        for j, res2 in enumerate(traj.topology.residues):
            if i == j: continue
            atom2_coords = residue_com(traj, res2)
            dist = np.linalg.norm(atom2_coords - atom1_coords)
            if dist <= local_dist_cutoff:
                close_atom_indices.append(j)
        
        traj.superpose(traj, atom_indices=close_atom_indices, ref_atom_indices=close_atom_indices)
        #positions = traj.xyz[:,i,:] - traj.xyz[0,i,:]
        atom1_coords_aligned = residue_com(traj, res1)
        positions = np.zeros((traj.n_frames, 3))
        for L in range(traj.n_frames):
            positions[L,:] = residue_com(traj, res1, frame=L) - atom1_coords_aligned
        
        locally_aligned_nodes[i, 0:num_frames] = positions[:,0]
        locally_aligned_nodes[i, num_frames:2*num_frames] = positions[:,1]
        locally_aligned_nodes[i, 2*num_frames:3*num_frames] = positions[:,2]
    
    return locally_aligned_nodes

def run_netcalc(graph, locally_aligned_nodes, num_nodes, num_node_pairs, silent=False):
    # Output correlation matrix
    R = cuarray.FloatCuArray()

    # Define all pair correlations that will be computed
    # Correlation pairs to compute
    ab = cuarray.IntCuArray()
    ab.init(num_node_pairs, 2)
    for i in range(num_nodes):
        for j in range(num_nodes):
            node_pair_index = i*num_nodes + j
            ab[node_pair_index][0] = i
            ab[node_pair_index][1] = j
            
    # Number of data points
    n = graph.numFrames()

    # Dimensionality of the data
    d = 3


    xd = 2

    # K-nearest-neighbors
    k = 6

    # CUDA platform
    platform = 0
        
    graph.nodeCoordinates().fromNumpy2D(locally_aligned_nodes.astype(np.float32))

    # Compute generalized correlation and output to proteinG_R
    if not silent:
        print("Performing generalized correlation computation "\
          f"on {n} data points with {num_nodes} nodes.")
    netcalc.generalizedCorrelation(
        X=graph.nodeCoordinates(),
        R=R,
        ab=ab,
        k=k,
        n=n,
        d=d,
        xd=xd,
        platform=platform,
    )

    # Gen. Corr. in numpy array object
    R_np = R.toNumpy2D().reshape(num_nodes, num_nodes)
    return R_np
    
def write_matrix(R_np, corr_matrix_filename, silent=False):
    if not silent:
        print("Saving matrix to:", corr_matrix_filename)
    np.savetxt(corr_matrix_filename, R_np)
    return

if __name__ == "__main__":
    trajectory_file = "1pgb_stripped.dcd"
    topology_file = "1pgb_stripped.pdb"
    starttime = time.time()
    traj, first_frame, last_frame, num_frames = load_trajectory(
        topology_file, trajectory_file)
    graph, num_nodes, num_node_pairs = create_netchem_graph(
        topology_file, trajectory_file, first_frame, last_frame)
    locally_aligned_nodes = locally_align_trajectory(
        traj, first_frame, last_frame, num_frames, num_nodes)
    R_np = run_netcalc(graph, locally_aligned_nodes, num_nodes, num_node_pairs)
    corr_matrix_filename = "proteinG_corr_matrix.txt"
    write_matrix(R_np, corr_matrix_filename)
    print(f"Total time: {time.time()-starttime:.3f} s.")
