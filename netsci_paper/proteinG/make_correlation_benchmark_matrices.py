"""
Find the 2D generalized correlation matrices for the B1 domain
of protein G multiple times, and use the JackKnife+ method to
estimate bias and noise within the dataset.
"""
import time
import os

import numpy as np

import make_correlation_matrix

NUM_MODELS = 100
TEMP_DCD = "DELETE_ME.dcd"


if __name__ == "__main__":
    trajectory_file = "1pgb_stripped.dcd"
    topology_file = "1pgb_stripped.pdb"
    starttime = time.time()
    traj, first_frame, last_frame, num_frames = make_correlation_matrix.load_trajectory(
        topology_file, trajectory_file, silent=True)
    
    for i in range(NUM_MODELS):
        print("Model:", i)
        indices = np.random.choice(np.arange(num_frames), size=num_frames, replace=True)
        subtraj = traj[indices]
        subtraj.save_dcd(TEMP_DCD)
        graph, num_nodes, num_node_pairs = make_correlation_matrix.create_netchem_graph(
            topology_file, TEMP_DCD, first_frame, last_frame, silent=True)
        locally_aligned_nodes = make_correlation_matrix.locally_align_trajectory(
            subtraj, first_frame, last_frame, num_frames, num_nodes, silent=True)
        R_np = make_correlation_matrix.run_netcalc(graph, locally_aligned_nodes, num_nodes, num_node_pairs, silent=True)
        corr_matrix_filename = f"proteinG_corr_bootstrap_{i}.txt"
        make_correlation_matrix.write_matrix(R_np, corr_matrix_filename)
    
    print(f"Total time: {time.time()-starttime:.3f} s.")
    os.remove(TEMP_DCD)
