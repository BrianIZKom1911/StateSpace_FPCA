#%%
import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval, legvander

freq = 'month' # 'month' or 'year' # change here for different frequency and folder
thisdir = os.path.dirname(__file__)
savedir = os.path.abspath(os.path.join(thisdir, '..', 'output', freq))
b_array = np.loadtxt(os.path.join(savedir, f'list_ridge1betas100_{freq}_ystd.csv'), delimiter=',')
#%%
# Fit PCA
pca = PCA()
Z_scores = pca.fit_transform(b_array) # (275x100) # each column is a PC time series
explained_var = pca.explained_variance_ratio_ # length 100
components = pca.components_ # (100x100) loading vectors w

np.savetxt(os.path.join(savedir, 'Z_scores_cov_100.csv'), Z_scores, delimiter=',')
np.savetxt(os.path.join(savedir, 'weights_cov_100.csv'), components, delimiter=',')
np.savetxt(os.path.join(savedir, 'explained_var_100.csv'), explained_var, delimiter=',')
#%%
# Scree plot - explained variance
plt.plot(np.cumsum(explained_var))
plt.xlabel("Number of PCs")
plt.ylabel("Cumulative explained variance")
plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(savedir, 'fig_var_explain.png'), dpi=300, bbox_inches='tight')
#%%
# Find the linear combination of basis functions for leading PCs (5)
# Simulate temperature curves
x = np.linspace(-5, 35, 80)
x_grid = 2*(x+10)/50 - 1 # map to [-1, 1] the same way as real data
X_poly = legvander(x_grid, 99) # (80x100)
q = 5 # number of leading PCs
Wq = components[:q, :].T # (100x5)
T = len(b_array)
# Faster way: Use broadcasting to compute all at once
Psi = X_poly @ Wq

# Save the data frame
col_names = ['r']+[f'PC{i+1}' for i in np.arange(Psi.shape[1])]
df = pd.DataFrame(np.hstack((x.reshape(-1, 1), Psi)), columns=col_names)
df.to_csv(os.path.join(savedir, 'df_leadingPCs5.csv'), index=None)
# End of script