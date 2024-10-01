"""
Produce a plot that resembles Fig. 2 of the Kraskov 2008 paper.
"""

import numpy as np
import matplotlib.pyplot as plt

cov_values = [0.9, 0.6, 0.3, 0.0] #results_matrix[1:,0]
colors = ["k", "r", "g", "b"]

SLICE = 11

def read_list_file(filename, bootstrap=False):
    inv_n_values_list = []
    average_list = []
    std_list = []
    on_index = -1
    with open(filename, "r") as f:
        for line in f.readlines():
            if line.startswith("# N"):
                N = float(line.strip().split()[2])
                inv_n_values_list.append(1 / N)
                if on_index > -1:
                    if bootstrap:
                        mean = np.mean(data_list) #np.median(data_list)
                        average_list.append(mean)
                        diff_list = []
                        for value in data_list:
                            diff = value - mean
                            diff_list.append(diff)

                        error = np.std(diff_list)
                        std_list.append(1.0*error)
                    else:
                        average_list.append(np.mean(data_list))
                        std_list.append(1.0*np.std(data_list))
                        
                on_index += 1
                data_list = []
                
            elif line.startswith("#"):
                continue
                
            else:
                data_list.append(float(line.strip()))
    
    if bootstrap:
        mean = np.mean(data_list) #np.median(data_list)
        average_list.append(mean)
        diff_list = []
        for value in data_list:
            diff = value - mean
            diff_list.append(diff)

        error = np.std(diff_list)
        std_list.append(1.0*error)
    else:
        average_list.append(np.mean(data_list))
        std_list.append(1.0*np.std(data_list))
        
    return inv_n_values_list, average_list, std_list

cov0_list_filename = "gaussian_I2_minus_Iexact0.0_list.txt"
cov3_list_filename = "gaussian_I2_minus_Iexact0.3_list.txt"
cov6_list_filename = "gaussian_I2_minus_Iexact0.6_list.txt"
cov9_list_filename = "gaussian_I2_minus_Iexact0.9_list.txt"

cov0_boot_filename = "gaussian_I2_minus_Iexact0.0_bootstrap.txt"
cov3_boot_filename = "gaussian_I2_minus_Iexact0.3_bootstrap.txt"
cov6_boot_filename = "gaussian_I2_minus_Iexact0.6_bootstrap.txt"
cov9_boot_filename = "gaussian_I2_minus_Iexact0.9_bootstrap.txt"

inv_n_values_list0, average_list0, std_list0 = read_list_file(
    cov0_list_filename)
inv_n_values_list3, average_list3, std_list3 = read_list_file(
    cov3_list_filename)
inv_n_values_list6, average_list6, std_list6 = read_list_file(
    cov6_list_filename)
inv_n_values_list9, average_list9, std_list9 = read_list_file(
    cov9_list_filename)
    
inv_n_values_boot0, average_boot0, std_boot0 = read_list_file(
    cov0_boot_filename, bootstrap=True)
inv_n_values_boot3, average_boot3, std_boot3 = read_list_file(
    cov3_boot_filename, bootstrap=True)
inv_n_values_boot6, average_boot6, std_boot6 = read_list_file(
    cov6_boot_filename, bootstrap=True)
inv_n_values_boot9, average_boot9, std_boot9 = read_list_file(
    cov9_boot_filename, bootstrap=True)

list_zeros0 = np.zeros(len(inv_n_values_list0))
list_zeros3 = np.zeros(len(inv_n_values_list3))
list_zeros6 = np.zeros(len(inv_n_values_list6))
list_zeros9 = np.zeros(len(inv_n_values_list9))
boot_zeros0 = np.zeros(len(inv_n_values_boot0))
boot_zeros3 = np.zeros(len(inv_n_values_boot3))
boot_zeros6 = np.zeros(len(inv_n_values_boot6))
boot_zeros9 = np.zeros(len(inv_n_values_boot9))

central_limit = np.zeros(len(inv_n_values_list0))
for i, inv_N in enumerate(inv_n_values_list0):
    central_limit[i] = np.sqrt(inv_N)

fig, axs = plt.subplots(4, 2)
axs[0][0].errorbar(inv_n_values_list0, list_zeros0, std_list0, fmt="k.-", label=f"r = 0.0", capsize=2)
axs[1][0].errorbar(inv_n_values_list3, list_zeros3, std_list3, fmt="k.-", label=f"r = 0.3", capsize=2)
axs[2][0].errorbar(inv_n_values_list6, list_zeros6, std_list6, fmt="k.-", label=f"r = 0.6", capsize=2)
axs[3][0].errorbar(inv_n_values_list9, list_zeros9, std_list9, fmt="k.-", label=f"r = 0.9", capsize=2)

axs[0][0].errorbar(inv_n_values_boot0, boot_zeros0, std_boot0, fmt="r.-", label=f"r = 0.0", capsize=2)
axs[1][0].errorbar(inv_n_values_boot3, boot_zeros3, std_boot3, fmt="r.-", label=f"r = 0.3", capsize=2)
axs[2][0].errorbar(inv_n_values_boot6, boot_zeros6, std_boot6, fmt="r.-", label=f"r = 0.6", capsize=2)
axs[3][0].errorbar(inv_n_values_boot9, boot_zeros9, std_boot9, fmt="r.-", label=f"r = 0.9", capsize=2)

axs[0][0].errorbar(inv_n_values_boot0, central_limit, np.zeros(central_limit.shape), fmt="g-")
axs[1][0].errorbar(inv_n_values_boot0, central_limit, np.zeros(central_limit.shape), fmt="g-")
axs[2][0].errorbar(inv_n_values_boot0, central_limit, np.zeros(central_limit.shape), fmt="g-")
axs[3][0].errorbar(inv_n_values_boot0, central_limit, np.zeros(central_limit.shape), fmt="g-")

axs[0][1].errorbar(inv_n_values_list0[SLICE:], list_zeros0[SLICE:], std_list0[SLICE:], fmt="k.-", label=f"r = 0.0", capsize=2)
axs[1][1].errorbar(inv_n_values_list3[SLICE:], list_zeros3[SLICE:], std_list3[SLICE:], fmt="k.-", label=f"r = 0.3", capsize=2)
axs[2][1].errorbar(inv_n_values_list6[SLICE:], list_zeros6[SLICE:], std_list6[SLICE:], fmt="k.-", label=f"r = 0.6", capsize=2)
axs[3][1].errorbar(inv_n_values_list9[SLICE:], list_zeros9[SLICE:], std_list9[SLICE:], fmt="k.-", label=f"r = 0.9", capsize=2)

axs[0][1].errorbar(inv_n_values_boot0[SLICE:], boot_zeros0[SLICE:], std_boot0[SLICE:], fmt="r.-", label=f"r = 0.0", capsize=2)
axs[1][1].errorbar(inv_n_values_boot3[SLICE:], boot_zeros3[SLICE:], std_boot3[SLICE:], fmt="r.-", label=f"r = 0.3", capsize=2)
axs[2][1].errorbar(inv_n_values_boot6[SLICE:], boot_zeros6[SLICE:], std_boot6[SLICE:], fmt="r.-", label=f"r = 0.6", capsize=2)
axs[3][1].errorbar(inv_n_values_boot9[SLICE:], boot_zeros9[SLICE:], std_boot9[SLICE:], fmt="r.-", label=f"r = 0.9", capsize=2)

axs[0][1].errorbar(inv_n_values_boot0[SLICE:], central_limit[SLICE:], np.zeros(central_limit[SLICE:].shape), fmt="g-")
axs[1][1].errorbar(inv_n_values_boot0[SLICE:], central_limit[SLICE:], np.zeros(central_limit[SLICE:].shape), fmt="g-")
axs[2][1].errorbar(inv_n_values_boot0[SLICE:], central_limit[SLICE:], np.zeros(central_limit[SLICE:].shape), fmt="g-")
axs[3][1].errorbar(inv_n_values_boot0[SLICE:], central_limit[SLICE:], np.zeros(central_limit[SLICE:].shape), fmt="g-")

fig.supylabel("$I^{(2)}(X,Y) + \\frac{1}{2} log(1-r^2)$")
axs[0][0].set_ylabel("r = 0.0")
axs[1][0].set_ylabel("r = 0.3")
axs[2][0].set_ylabel("r = 0.6")
axs[3][0].set_ylabel("r = 0.9")
axs[3][0].set_xlabel("1/N")
axs[3][1].set_xlabel("1/N")
#fig.title("Gaussian System $I^{(2)}(X,Y)-I_{exact}(X,Y)$")
#plt.ylim([0, 0.05])
#plt.ylim([-0.015, 0.01])
plt.tick_params(direction="in", right="on", top="on")
#plt.legend()
plt.tight_layout()

plt.show()
