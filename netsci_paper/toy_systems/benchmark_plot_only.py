"""
Benchmark the GPU vs. CPU implementations.
"""

import time

import numpy as np

import matplotlib.pyplot as plt

"""
def plot_GPU_vs_CPU_times(N_values, GPU_times_per_N, CPU_times_per_N, filename=None):
    pi_fig, ax1 = plt.subplots()
    ax1.plot(N_values, GPU_times_per_N, "r", marker="o", linestyle="-",
             label="GPU benchmark")
    ax1.plot(N_values, CPU_times_per_N, "g", marker="o", linestyle="-",
             label="CPU benchmark")
    ax1.set_ylabel("Time per MI calculation (s)")
    ax1.set_xlabel("Number of realizations")
    
    plt.yscale("log")
    plt.xscale("log")
    plt.title("GPU and CPU benchmarks")
    plt.legend()
    plt.tight_layout()
    
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
"""

def plot_GPU_vs_CPU_times(N_values, 
        GPU1_times_per_N_10, GPU2_times_per_N_10, CPU_times_per_N_10, 
        GPU1_times_per_N_30, GPU2_times_per_N_30, CPU_times_per_N_30, 
        GPU1_times_per_N_100, GPU2_times_per_N_100, CPU_times_per_N_100, 
        filename=None):
    pi_fig, ax1 = plt.subplots()
    ax1.plot(N_values, GPU1_times_per_N_10, "lightcoral", marker="o", linestyle="-",
             label="GTX1080 10 nodes")
    ax1.plot(N_values, GPU1_times_per_N_30, "indianred", marker="o", linestyle="-",
             label="GTX1080 30 nodes")
    ax1.plot(N_values, GPU1_times_per_N_100, "brown", marker="o", linestyle="-",
             label="GTX1080 100 nodes")
    ax1.plot(N_values, GPU2_times_per_N_10, "limegreen", marker="o", linestyle="-",
             label="RTX6000 Ada 10 nodes")
    ax1.plot(N_values, GPU2_times_per_N_30, "forestgreen", marker="o", linestyle="-",
             label="RTX6000 Ada 30 nodes")
    ax1.plot(N_values, GPU2_times_per_N_100, "darkgreen", marker="o", linestyle="-",
             label="RTX6000 Ada 100 nodes")
    ax1.plot(N_values, CPU_times_per_N_10, "lightsteelblue", marker="o", linestyle="-",
             label="CPU 10 nodes")
    ax1.plot(N_values, CPU_times_per_N_30, "cornflowerblue", marker="o", linestyle="-",
             label="CPU 30 nodes")
    ax1.plot(N_values, CPU_times_per_N_100, "royalblue", marker="o", linestyle="-",
             label="CPU 100 nodes")
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

if __name__ == "__main__":
    
    starttime = time.time()
    N_values = [31, 100, 316, 1000, 3160, 10000, 31600]
    #GPU_times_per_N, CPU_times_per_N = benchmark_GPU_vs_CPU_2x2(N_values)
    #GPU_times_per_N = [ 0.000542, 0.000179, 6.333e-07, 2.672e-06, 1.262e-05, 4.053e-05]
    GPU1_times_per_N_10 = np.loadtxt("GPU_gtx1080_times_M_10.txt")
    GPU1_times_per_N_30 = np.loadtxt("GPU_gtx1080_times_M_30.txt")
    GPU1_times_per_N_100 = np.loadtxt("GPU_gtx1080_times_M_100.txt")
    
    GPU2_times_per_N_10 = np.loadtxt("GPU_rtx6000ada_times_M_10.txt")
    GPU2_times_per_N_30 = np.loadtxt("GPU_rtx6000ada_times_M_30.txt")
    GPU2_times_per_N_100 = np.loadtxt("GPU_rtx6000ada_times_M_100.txt")
    
    CPU_times_per_N_10 = np.loadtxt("CPU_times_M_10.txt")
    CPU_times_per_N_30 = np.loadtxt("CPU_times_M_30.txt")
    CPU_times_per_N_100 = np.loadtxt("CPU_times_M_100.txt")
    
    #CPU_times_per_N_10 = [6.641e-05, 0.000232, 0.000817,  0.00292, 0.0101, 0.0363, np.nan]
    #CPU_times_per_N_30 = [ 0.000881, 0.00320, 0.0114,  0.0406,  0.143, 0.502, np.nan]
    #CPU_times_per_N_100 = [0.00681, 0.0244, 0.088, 0.323, 1.153, 4.073, np.nan]
    
    plot_GPU_vs_CPU_times(N_values, 
        GPU1_times_per_N_10, GPU2_times_per_N_10, CPU_times_per_N_10, 
        GPU1_times_per_N_30, GPU2_times_per_N_30, CPU_times_per_N_30, 
        GPU1_times_per_N_100, GPU2_times_per_N_100, CPU_times_per_N_100)
    
    #CPU_times_per_N = [6.641e-05, 0.000232, 0.000817,  0.00292, 0.0101, 0.0363]
    
    """
    plot_GPU_vs_CPU_times(N_values, GPU_times_per_N, CPU_times_per_N, filename="benchmark_2x2.png")
    print("Time to do M=10 (s):", time.time() - starttime)
    
    #GPU_times_per_N, CPU_times_per_N = benchmark_GPU_vs_CPU_MxM(M, N_values)
    GPU_times_per_N = [ 0.000297,  9.047e-05, 3.666e-05, 3.020e-05, 0.000133, 0.000600, 0.00203]
    CPU_times_per_N = [ 0.000881, 0.00320, 0.0114,  0.0406,  0.143, 0.502, np.nan]
    #plot_GPU_vs_CPU_times(N_values, GPU_times_per_N, CPU_times_per_N, filename="benchmark_10x10.png")
    print("Time to do M=30 (s):", time.time() - starttime)
    
    N_values = [31, 100, 316, 1000, 3160, 10000, 31600]
    M = 30
    #GPU_times_per_N, CPU_times_per_N = benchmark_GPU_vs_CPU_MxM(M, N_values)
    GPU_times_per_N = [0.189, 0.055, 0.0186, 0.00725, 0.00347, 0.00275, 0.0130, 0.0603]
    CPU_times_per_N = [0.00681, 0.0244, 0.088, 0.323, 1.153, 4.073, np.nan, np.nan]
    plot_GPU_vs_CPU_times(N_values, GPU_times_per_N, CPU_times_per_N, filename="benchmark_100x100.png")
    print("Time to do M=100 (s):", time.time() - starttime)
    """
