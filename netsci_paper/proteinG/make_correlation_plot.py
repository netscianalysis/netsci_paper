"""
Make the plot of the proteinG correlation matrix.
"""
import numpy as np
import matplotlib.pyplot as plt

R_np = np.loadtxt("proteinG_corr_matrix_local2.txt")
R_np1 = np.loadtxt("proteinG_pearson_matrix.txt")

# Ranges of data
num_nodes = 56

for i in range(num_nodes):
    for j in range(i, num_nodes):
        R_np[i,j] = R_np1[i,j]


R_np = np.flip(R_np, axis=0)
R_np1 = np.flip(R_np1, axis=0)


R_figure_x = [i for i in range(num_nodes)]
R_figure_y = [i for i in range(num_nodes)]

im = plt.imshow(R_np, vmin=0.0, vmax=1.0, extent=[1, num_nodes+1, 1, num_nodes+1],
                cmap=plt.cm.jet)
im.set_interpolation('bilinear')
plt.xlabel("Residue number")
plt.ylabel("Residue number")
cbar = plt.colorbar(im)
plt.show()
