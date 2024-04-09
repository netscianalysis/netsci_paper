"""
Benchmark the GPU vs. CPU implementations.
"""
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

import netcalc
import cuarray
import gaussian_mi


def benchmark_MxM(platform, M, N):
    k = 1
    d = 1
    xd = 2
    n_realizations = 1
    
    ab = cuarray.IntCuArray()
    num_nodes = M
    num_node_pairs = num_nodes**2

    # Define all pair correlations that will be computed
    ab.init(num_node_pairs, 2)
    for i in range(num_nodes):
        for j in range(num_nodes):
            node_pair_index = i*num_nodes + j
            ab[node_pair_index][0] = i
            ab[node_pair_index][1] = j
    
    gaussian_2D = gaussian_mi.make_gaussian_independent_data(M, N)
    X = cuarray.FloatCuArray()
    X.fromNumpy2D(gaussian_2D)
    I = cuarray.FloatCuArray()
    starttime = time.time()
    netcalc.mutualInformation(X, I, ab, k, N, xd, d, platform)
    total_time = time.time() - starttime
        
    return total_time



def benchmark_GPU_vs_CPU_MxM(M, N_values, attempt_N_with_GPU, attempt_N_with_CPU, gpu_only):
    GPU_times_per_N = []
    CPU_times_per_N = []
    for i, N in enumerate(N_values):
        if attempt_N_with_GPU[i]:
            print(f"Running GPU benchmark. N = {N}")
            GPU_time = benchmark_MxM(netcalc.GPU_PLATFORM, M, N)
            GPU_time_per_N = GPU_time / (0.5 * N * M * (M-1))
            print("Time per N (s):", GPU_time_per_N)
            GPU_times_per_N.append(GPU_time_per_N)
        else:
            print(f"Skipping N = {N} for GPU")
            GPU_times_per_N.append(np.nan)
        
        if not gpu_only:
            if attempt_N_with_CPU[i]:
                print(f"Running CPU benchmark. N = {N}")
                CPU_time = benchmark_MxM(netcalc.CPU_PLATFORM, M, N)
                CPU_time_per_N = CPU_time / (0.5 * N * M * (M-1))
                print("Time per N (s):", CPU_time / N)
                CPU_times_per_N.append(CPU_time_per_N)
            else:
                print(f"Skipping N = {N} for CPU")
                CPU_times_per_N.append(np.nan)
        
    return GPU_times_per_N, CPU_times_per_N

def plot_GPU_vs_CPU_times(N_values, 
        GPU_times_per_N_10, CPU_times_per_N_10, 
        GPU_times_per_N_100, CPU_times_per_N_100, 
        GPU_times_per_N_1000, CPU_times_per_N_1000, 
        filename=None):
    pi_fig, ax1 = plt.subplots()
    ax1.plot(N_values, GPU_times_per_N_10, "lightcoral", marker="o", linestyle="-",
             label="GPU 10 nodes")
    ax1.plot(N_values, GPU_times_per_N_100, "indianred", marker="o", linestyle="-",
             label="GPU 100 nodes")
    ax1.plot(N_values, GPU_times_per_N_1000, "brown", marker="o", linestyle="-",
             label="GPU 1000 nodes")
    ax1.plot(N_values, CPU_times_per_N_10, "lightsteelblue", marker="o", linestyle="-",
             label="CPU 10 nodes")
    ax1.plot(N_values, CPU_times_per_N_100, "cornflowerblue", marker="o", linestyle="-",
             label="CPU 100 nodes")
    ax1.plot(N_values, CPU_times_per_N_1000, "royalblue", marker="o", linestyle="-",
             label="CPU 1000 nodes")
    ax1.set_ylabel("Time per MI calculation (s)")
    ax1.set_xlabel("Number of data points")
    
    plt.yscale("log")
    plt.xscale("log")
    plt.title("GPU and CPU benchmarks")
    plt.legend()
    plt.tight_layout()
    
    if filename is None:
        plt.show()
    else:
        print("saving figure:", filename)
        plt.savefig(filename)

def save_results(GPU_times, CPU_times, suffix, gpu_only):
    np.savetxt(f"GPU_times_{suffix}.txt", GPU_times)
    if not gpu_only:
        np.savetxt(f"CPU_times_{suffix}.txt", CPU_times)
    return
    

if __name__ == "__main__":
    
    if "gpu_only" in sys.argv:
        gpu_only = True
    else:
        gpu_only = False
    
    #N_values = [100, 316, 1000, 3160, 10000, 31600, 100000]
    N_values = [31, 100, 316, 1000, 3160, 10000, 31600]
    #N_values = [10, 32, 100, 316]
    starttime = time.time()
    M = 10
    attempt_N_with_GPU = [True, True, True, True, True, True, True]
    #attempt_N_with_GPU = [True, True, True, True]
    attempt_N_with_CPU = [True, True, True, True, True, False, False]
    #attempt_N_with_CPU = [True, True, True, False]
    GPU_times_per_N_10, CPU_times_per_N_10 = benchmark_GPU_vs_CPU_MxM(M, N_values, attempt_N_with_GPU, attempt_N_with_CPU, gpu_only)
    print("Time to do M=10 (s):", time.time() - starttime)
    save_results(GPU_times_per_N_10, CPU_times_per_N_10, "M_10", gpu_only)
    
    M = 30
    attempt_N_with_GPU = [True, True, True, True, True, True, True]
    #attempt_N_with_GPU = [True, True, True, True]
    attempt_N_with_CPU = [True, True, True, True, False, False, False]
    #attempt_N_with_CPU = [True, True, False, False]
    GPU_times_per_N_30, CPU_times_per_N_30 = benchmark_GPU_vs_CPU_MxM(M, N_values, attempt_N_with_GPU, attempt_N_with_CPU, gpu_only)
    print("Time to do M=30 (s):", time.time() - starttime)
    save_results(GPU_times_per_N_30, CPU_times_per_N_30, "M_30", gpu_only)
    
    M = 100
    attempt_N_with_GPU = [True, True, True, True, True, True, True]
    #attempt_N_with_GPU = [True, True, True, False]
    attempt_N_with_CPU = [True, True, True, False, False, False, False]
    #attempt_N_with_CPU = [True, True, False, False]
    GPU_times_per_N_100, CPU_times_per_N_100 = benchmark_GPU_vs_CPU_MxM(M, N_values, attempt_N_with_GPU, attempt_N_with_CPU, gpu_only)
    print("Time to do M=100 (s):", time.time() - starttime)
    save_results(GPU_times_per_N_100, CPU_times_per_N_100, "M_100", gpu_only)
    
    #plot_GPU_vs_CPU_times(N_values, GPU_times_per_N_10, CPU_times_per_N_10, 
    #                      GPU_times_per_N_100, CPU_times_per_N_100, 
    #                      GPU_times_per_N_1000, CPU_times_per_N_1000, filename="benchmark_all.png")
