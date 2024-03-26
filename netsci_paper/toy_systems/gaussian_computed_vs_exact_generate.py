"""
Recreate the data to generate a replica of Fig. 2 from Kraskov 2008
origin MI paper.
"""

import time

import netcalc
import gaussian_mi

if __name__ == "__main__":
    k = 1
    d = 1
    xd = 2
    platform = netcalc.GPU_PLATFORM
    covariance_values = [0.9, 0.6, 0.3, 0.0]
    N_values = [20, 25, 30, 35, 40, 45, 50, 60, 80, 100, 140, 180, 220, 280,
                400, 800, 1600, 6400, 20000, 40000]
    starttime = time.time()
    benchmark = gaussian_mi.run_gaussian_mi_all(covariance_values, N_values,
        k, xd, d, platform, filename_prefix="gaussian_I2_minus_Iexact")
    print(f"Total time: {time.time()-starttime:.3f}")
