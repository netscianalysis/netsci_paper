"""
Make the plot of the proteinG correlation matrix.
"""
import numpy as np
import matplotlib.pyplot as plt

# Ranges of data
num_nodes = 56

#R_median = np.zeros((num_nodes, num_nodes))

R_matrix_list = []

for k in range(100):
    filename = f"proteinG_corr_bootstrap_{k}.txt"
    R_np_i = np.loadtxt(filename)
    R_np_i = np.flip(R_np_i, axis=0)
    R_matrix_list.append(R_np_i)

R_mean = np.mean(R_matrix_list, axis=0)

#im = plt.imshow(R_mean, vmin=0.0, extent=[1, num_nodes+1, 1, num_nodes+1],
#                cmap=plt.cm.jet)
#im.set_interpolation('bilinear')
#plt.xlabel("Residue number")
#plt.ylabel("Residue number")
#cbar = plt.colorbar(im)
#plt.show()

R_q_list = []
for R_np_i in R_matrix_list:
    R_q = R_np_i - R_mean
    R_q_list.append(R_q)

R_error = np.std(R_q_list, axis=0)

R_error_10 = np.std(R_q_list[:10], axis=0)

for i in range(num_nodes):
    for j in range(num_nodes - i):
        R_error[i,j] = R_error_10[i,j]

print("R_error.shape:", R_error.shape)

#cmap = plt.cm.cool
#cmap = plt.cm.PiYG
cmap = plt.cm.RdPu
#cmap = plt.cm.jet
vmax = None
#vmax = 1.0

im = plt.imshow(R_error, vmin=0.0, vmax=vmax, extent=[1, num_nodes+1, 1, num_nodes+1],
                cmap=cmap)
im.set_interpolation('bilinear')
plt.xlabel("Residue number")
plt.ylabel("Residue number")
cbar = plt.colorbar(im)
plt.show()
