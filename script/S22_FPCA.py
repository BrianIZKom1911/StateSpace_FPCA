#%%
import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval, legvander
#%% Define procedure
def DOTHIS_PCA(q1, q2, mx_obs, dir_save, savefile=False):
    # Fit PCA
    pca = PCA(n_components=q1) # e.g. q1=20, q2=5
    Z_scores = pca.fit_transform(mx_obs) # (275x20) # each column is a PC time series
    explained_var = pca.explained_variance_ratio_ # length 20
    components = pca.components_ # (20x100) loading vectors w

    # Find the linear combination of basis functions for leading PCs (5)
    # Simulate temperature curves
    x = np.linspace(-5, 35, 80)
    x_grid = 2*(x+10)/50 - 1 # map to [-1, 1] the same way as real data
    X_poly = legvander(x_grid, 99) # (80x100)
    Aq = components[:q2, :].T # (5x100) --> (100x5)
    T = len(mx_obs)
    Psi = X_poly @ Aq # use broadcasting to compute all at once

    if savefile:
        # Save PCA results
        np.savetxt(os.path.join(dir_save, f'Z_scores_cov_{q1}.csv'), Z_scores, delimiter=',')
        np.savetxt(os.path.join(dir_save, f'weights_cov_{q1}.csv'), components, delimiter=',')
        np.savetxt(os.path.join(dir_save, f'explained_var_{q1}.csv'), explained_var, delimiter=',')        
        # Save the dataframe of FPCs
        col_names = ['r']+[f'PC{i+1}' for i in np.arange(Psi.shape[1])]
        df = pd.DataFrame(np.hstack((x.reshape(-1, 1), Psi)), columns=col_names)
        df.to_csv(os.path.join(dir_save, f'df_leadingPCs{q2}.csv'), index=None)

#%%
thisdir = os.path.dirname(__file__)
outputdir = os.path.abspath(os.path.join(thisdir, '..', 'output'))
subfolders = ['1.0_month', '1.1_month', '1.0_year', '1.1_year', '2_year']
#df_std = pd.read_csv(os.path.join(outputdir, 'df_std.csv'), index_col=False)
#df_std['key'] = subfolders

inputfiles = [f'{s}_ridgecoefs100' for s in subfolders]
for i, subfd in enumerate(subfolders):
    infe = inputfiles[i]
    savedir2 = os.path.join(outputdir, subfd)
    os.makedirs(savedir2, exist_ok=True)

    matches = [f for f in os.listdir(outputdir) if infe in f]
    inpath = os.path.join(outputdir, matches[0])
    c_array = np.loadtxt(inpath, delimiter=',')

    #mu = df_std.loc[df_std['key']==subfd, 'mu'].iat[0]
    #sd = df_std.loc[df_std['key']==subfd, 'sd'].iat[0]
    # DO THIS
    DOTHIS_PCA(20, 5, c_array, savedir2, savefile=True)

# End of script.