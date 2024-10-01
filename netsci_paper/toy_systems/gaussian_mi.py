"""
Base functions for computing mutual information (MI) of a 2D Gaussian
distribution.
"""

import time
import numpy as np
import cuarray
import netchem
import netcalc

# TODO: move to a common module for Gaussian calcs
def make_covariance_matrix(covariance_value):
    cov = np.array([[1.0, covariance_value], [covariance_value, 1.0]])
    return cov

# TODO: move to a common module for Gaussian calcs
def compute_I_exact(covariance_value):
    cov = make_covariance_matrix(covariance_value)
    mi_control = -0.5*np.log(np.linalg.det(cov))
    return mi_control

# TODO: move to a common module for Gaussian calcs
def make_gaussian_2D_points(N, covariance_value):
    gaussian_2D_mean = np.zeros(2)
    gaussian_2D_cov = make_covariance_matrix(covariance_value)
    gaussian_2D = np.random.multivariate_normal(
        mean=gaussian_2D_mean,
        cov=gaussian_2D_cov,
        size=N,
    ).T.astype(np.float32)
    return gaussian_2D

def make_ab_for_mi():
    ab = cuarray.IntCuArray()
    ab.init(2, 2)
    ab[0][0] = 0
    ab[0][1] = 1
    ab[1][0] = 1
    ab[1][1] = 0
    return ab

def make_gaussian_independent_data(M, N):
    gaussian_2D = np.zeros((M, N)).astype(np.float32)
    for i in range(M):
        gaussian_2D[i,:] = np.random.normal(size=N)
    
    return gaussian_2D

def run_gaussian_mi(n_realizations, N, covariance_value, ab, k, xd, d, platform):
    #I_values_sum = 0.0
    I_values_list = []
    total_time = 0.0
    for kay in range(n_realizations):
        gaussian_2D = make_gaussian_2D_points(N, covariance_value)
        X = cuarray.FloatCuArray()
        X.fromNumpy2D(gaussian_2D)
        I = cuarray.FloatCuArray()
        starttime = time.time()
        netcalc.mutualInformation(X, I, ab, k, N, xd, d, platform)
        total_time += time.time() - starttime
        print("I[0][0]:", I[0][0])
        I_value = I[0][0]
        #assert np.isfinite(I_value), "NaN detected!"
        #I_values_sum += I_value
        I_values_list.append(I_value)
    
    return I_values_list, total_time
    #I_values_avg = I_values_sum / n_realizations
    #return I_values_avg, total_time

def run_gaussian_mi_all(covariance_values, N_values,
        k, xd, d, platform, filename_prefix="gaussian_I2_minus_Iexact"):
    ab = make_ab_for_mi()
    starttime = time.time()
    for i, covariance_value in enumerate(covariance_values):
        results_matrix = np.zeros((2, len(N_values)+1))
        print(f"covariance value: {covariance_value}, "
              f"{i+1} of {len(covariance_values)}, "
              f"time: {time.time()-starttime:.3f}")
        I_exact = compute_I_exact(covariance_value)
        #results_matrix[1, 0] = covariance_value
        list_filename = f"{filename_prefix}{covariance_value}_list.txt"
        with open(list_filename, "w") as f:
            f.write(f"# covariance_value: {covariance_value}\n")
            
        for j, N in enumerate(N_values):
            print(f"N value: {N}, {j+1} of len(N_values)," 
                  f"time: {time.time()-starttime:.3f}")
            if N <= 100:
                n_realizations = 200 #2000000
            elif N <= 1000:
                n_realizations = 200 #500000
            else:
                n_realizations = 200 #100000
                
            I_values_list, total_time = run_gaussian_mi(n_realizations, N, covariance_value, ab, k, xd, d, platform)
            I_values_avg = np.mean(I_values_list)
            results_matrix[1,j+1] = I_values_avg - I_exact
            
            with open(list_filename, "a") as f:
                f.write(f"# N: {N}\n")
                for I_value in I_values_list:
                    I_value_minus_I_exact = I_value - I_exact
                    f.write(f"{I_value_minus_I_exact}\n")
        
        for j, N in enumerate(N_values):
            results_matrix[0, j+1] = N
    
        np.savetxt(f"{filename_prefix}{covariance_value}.txt", results_matrix)
    
    return
