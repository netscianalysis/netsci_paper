"""
Make the 2D generalized correlation matrices for the B1 domain
of protein G, showing correlation between residues.
"""
import time

import numpy as np

import cuarray
import netchem
import netcalc

def pearson_correlation(
    node_coordinates: cuarray.FloatCuArray,
    num_nodes: int,
    num_frames: int,
) -> np.ndarray:
    node_coordinates_copy = cuarray.FloatCuArray()
    node_coordinates_copy.fromCuArray(
        cuArray=node_coordinates,
        start=0,
        end=num_nodes-1,
        m=node_coordinates.m(),
        n=node_coordinates.n(),
    )
    for i in range(num_nodes):
        u = cuarray.FloatCuArray()
        a = cuarray.FloatCuArray()
        a.fromCuArray(
            cuArray=node_coordinates_copy,
            start=i,
            end=i,
            m=3,
            n=num_frames
        )
        netcalc.mean(
            a=a,
            u=u,
            m=3,
            n=num_frames,
            platform=netcalc.GPU_PLATFORM,
        )
        for j in range(num_frames):
            node_coordinates_copy[i][j] -= u[0][0]
            node_coordinates_copy[i][j] -= u[0][1]
            node_coordinates_copy[i][j] -= u[0][2]
    np_node_coordinates = node_coordinates_copy.toNumpy2D().reshape(
        num_nodes,
        3,
        num_frames
    )
    X = np.zeros((num_nodes, num_frames))
    for i in range(num_nodes):
        for j in range(num_frames):
            X[i, j] = np.linalg.norm(np_node_coordinates[i, :, j])
    C = np.corrcoef(X)
    return C

starttime = time.time()
trajectory_file = "1pgb_stripped.dcd"
topology_file = "1pgb_stripped.pdb"
first_frame = 0
last_frame = 19499

print("Creating netchem graph")
graph = netchem.Graph()
graph.init(
    trajectoryFile=trajectory_file,
    topologyFile=topology_file,
    firstFrame=first_frame,
    lastFrame=last_frame,
)

# Output correlation matrix
R = cuarray.FloatCuArray()

# Correlation pairs to compute
ab = cuarray.IntCuArray()
num_nodes = graph.numNodes()
num_node_pairs = num_nodes**2

# Define all pair correlations that will be computed
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

# Compute generalized correlation and output to proteinG_R
print("Performing Pearson correlation computation "\
      f"on {n} data points with {num_nodes} nodes.")
print(f"Time: {time.time()-starttime:.3f} s.")


a = graph.nodeCoordinates() #.toNumpy2D()
#a = a.reshape(num_nodes, 3, last_frame+1)
r = pearson_correlation(a, num_nodes, 19500)

R_np = r.reshape(num_nodes, num_nodes)
pearson_matrix_filename = "proteinG_pearson_matrix.txt"
print("Saving matrix to:", pearson_matrix_filename)
np.savetxt(pearson_matrix_filename, R_np)
print(f"Total time: {time.time()-starttime:.3f} s.")
