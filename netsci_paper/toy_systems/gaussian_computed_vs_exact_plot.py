"""
Produce a plot that resembles Fig. 2 of the Kraskov 2008 paper.
"""

import numpy as np
import matplotlib.pyplot as plt

cov_values = [0.9, 0.6, 0.3, 0.0] #results_matrix[1:,0]
colors = ["k", "r", "g", "b"]

fig, ax = plt.subplots()
for i, cov_value in enumerate(cov_values):
    results_matrix = np.loadtxt(f"gaussian_I2_minus_Iexact{cov_value}.txt")
    N_values = results_matrix[0,1:]
    inv_N_values = N_values ** -1
    
    I_diff_values = results_matrix[1,1:]
    #fig.add_trace(go.Scatter(x=inv_N_values, y=I_diff_values, 
    #                         mode="lines+markers", name=f"r = {cov_value}"))
    ax.plot(inv_N_values, I_diff_values, f"{colors[i]}.-", label=f"r = {cov_value}")


plt.ylabel("$I^{(2)}(X,Y) + \\frac{1}{2} log(1-r^2)$")
plt.xlabel("1/N")
plt.title("Gaussian System $I^{(2)}(X,Y)-I_{exact}(X,Y)$")
plt.ylim([0, 0.05])
plt.ylim([-0.015, 0.01])
plt.tick_params(direction="in", right="on", top="on")
plt.legend()
plt.tight_layout()
plt.show()
