"""
Find a distribution of N Gaussians and use the Bootstrap 
method with replacement to estimate bias and noise within 
the dataset.
"""
import time
import os

import numpy as np
import cuarray
import netchem
import netcalc
import gaussian_mi

NUM_MODELS = 100
covariance_values = [0.9, 0.6, 0.3, 0.0]
N_values = [20, 25, 30, 35, 40, 45, 50, 60, 80, 100, 140, 180, 220, 280,
                400, 800, 1600, 6400, 20000, 40000]
k = 1
d = 1
xd = 2
platform = netcalc.GPU_PLATFORM
NUM_MODELS = 100
filename_prefix="gaussian_I2_minus_Iexact"

if __name__ == "__main__":
    
    ab = cuarray.IntCuArray()
    ab.init(2, 2)
    ab[0][0] = 0
    ab[0][1] = 1
    ab[1][0] = 1
    ab[1][1] = 0
    
    for i, covariance_value in enumerate(covariance_values):
        print("covariance_value:", covariance_value)
        list_filename = f"{filename_prefix}{covariance_value}_bootstrap.txt"
        with open(list_filename, "w") as f:
            f.write(f"# covariance_value: {covariance_value}\n")
        
        for j, N in enumerate(N_values):
            print("N:", N)
            gaussian_2D = gaussian_mi.make_gaussian_2D_points(N, covariance_value)
            I_values_list = []
            with open(list_filename, "a") as f:
                f.write(f"# N: {N}\n")
                
            for i in range(NUM_MODELS):
                frame_to_leave_out = np.random.choice(N)
                #print("frame_to_leave_out:", frame_to_leave_out)
                new_gaussian_2D = np.zeros(gaussian_2D.shape, dtype=np.float32)
                indices = np.random.choice(np.arange(N), size=N, replace=True)
                new_gaussian_2D[0,:] = gaussian_2D[0,indices]
                new_gaussian_2D[1,:] = gaussian_2D[1,indices]
                #new_gaussian_2D = np.delete(gaussian_2D, frame_to_leave_out, axis=1)
                X = cuarray.FloatCuArray()
                X.fromNumpy2D(new_gaussian_2D)
                I = cuarray.FloatCuArray()
                netcalc.mutualInformation(X, I, ab, k, N, xd, d, platform)
                #print("I[0][0]:", I[0][0])
                I_value = I[0][0]
                I_exact = gaussian_mi.compute_I_exact(covariance_value)
                I_value_minus_I_exact = I_value - I_exact
                with open(list_filename, "a") as f:
                    f.write(f"{I_value_minus_I_exact}\n")
            
