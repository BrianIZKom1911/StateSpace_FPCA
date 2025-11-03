# This file trims extremes inside the local method
# which is fine but not sophisticated as LOWESS
#%%
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

thisdir = os.path.dirname(__file__)
savedir = os.path.abspath(os.path.join(thisdir, '..', 'output', 'month'))
os.makedirs(savedir, exist_ok=True)
#%% # Import and process data
url = "https://raw.githubusercontent.com/BrianIZKom1911/electricload_CTRF/master/data_clean/NC_main.csv"
data = pd.read_csv(url, index_col=None)
# Adjust data type
data['datetime_UTC'] = pd.to_datetime(data['datetime_UTC'])
data['Year'] = data['Year'].astype(int)
data['Month'] = data['Month'].astype(int)
data['Day'] = data['Day'].astype(int)
data['Hour'] = data['Hour'].astype(int)
data['yyyymm'] = data['Year'].astype(str) + '-' + data['Month'].astype(str).str.zfill(2)
# Truncate extremes
dt1 = data[(data['temperature'] >= -10) & (data['temperature'] <= 40)]
#%% # I) Smooth the one-month load data using local linear kNN
def local_linear_knn(x, y, k, x_grid, lower_q=0.0, onesided=True, lambda_ridge=0.0):
    """
    Perform local linear regression using k-NN and Epanechnikov kernel weighting (09/02)
    
    Parameters:
    r (array-like): The independent variable
    y (array-like): The dependent variable
    k (int): Number of nearest neighbors to consider
    x_grid (array-like): Points at which to estimate the conditional expectation function
    lower_q (float in [0,1)): optional, lower quantile to drop within each neighborhood
    onesided (bool): optional, if True, also trim the upper quantile defined by 1-lower_q 
    lambda_ridge (float): optional, Ridge penalty

    Returns:
    np.ndarray: Estimated conditional expectation function at x_grid E[Y|X=x_grid]
    """
    x_norm = x.reshape(-1, 1)
    x_grid_norm = x_grid.reshape(-1, 1)
    
    # Find k-NN for each point in x_grid
    nbrs = NearestNeighbors(n_neighbors=k).fit(x_norm)
    distances, indices = nbrs.kneighbors(x_grid_norm)
    # Initialize output
    y_pred = np.zeros_like(x_grid)

    for i in range(len(x_grid)):
        # Find k-NN data points
        x_nbrs = x[indices[i]]
        y_nbrs = y[indices[i]]
        # Compute Epanechnikov weights
        dist_max = distances[i, -1] # The critical bandwidth for x_grid[i] (distance to its k-th nearest neighbor)
        u = distances[i] / dist_max # Normalized distances
        weights = 0.75*(1 - u**2) * (np.abs(u) <= 1) # Epanechnikov kernel weights

        # ---- trimming by neighborhood percentile ----
        if lower_q > 0:
            upper_q = 1 - lower_q
            ylb = np.quantile(y_nbrs, lower_q); yub = np.quantile(y_nbrs, upper_q)
            mask = y_nbrs >= ylb if onesided else (y_nbrs>=ylb) & (y_nbrs<=yub)
            x_nbrs = x_nbrs[mask]
            y_nbrs = y_nbrs[mask]
            weights = weights[mask]
        
        # Weighted least squares for local linear fit at x_grid[i]
        X = np.column_stack([np.ones_like(x_nbrs), x_nbrs - x_grid[i]])
        X_w = np.sqrt(weights)[:, None] * X
        y_w = np.sqrt(weights) * y_nbrs
        # Regularization
        beta = np.linalg.solve(
            X_w.T @ X_w + lambda_ridge * np.eye(X.shape[1]),
            X_w.T @ y_w
        )
        y_pred[i] = beta[0] # intercept = fitted value at x_grid[i]
    
    return y_pred

# Loop over year-month to average load for each period
months = dt1['yyyymm'].unique()
T = len(months) # 276 months
dt1['y_avg'] = np.nan

for t in range(T):
    mon = months[t]
    dt_m = dt1.loc[dt1['yyyymm'] == mon, ['yyyymm', 'load', 'temperature']]
    if len(dt_m) < 180: # skip months with too few data points
        continue
        
    # Local weighted linear rgr (kNN)
    r = dt_m['temperature'].values
    y = dt_m['load'].values
    kt = 10*np.sqrt(len(y)).astype(int) # Debatable (2)
    r0 = r.copy()
    y_avg = local_linear_knn(r, y, kt, r0, lower_q=0.015, onesided=False, lambda_ridge=0.1)
    
    # Store results
    dt1.loc[dt_m.index, 'y_avg'] = y_avg # assign output only to the rows that remain after filtering

# Save the average results
cols_save = ['datetime_UTC', 'Year', 'Month', 'Day', 'Hour', 'yyyymm', 'temperature', 'load', 'y_avg']
dt1[cols_save].to_csv(os.path.join(savedir, 'NC_avgload_month.csv'), index=False)
#%% # II) Fit polynomial basis to the monthly average load data
from sklearn.linear_model import Ridge
from numpy.polynomial.legendre import legval, legvander
from sklearn.linear_model import LinearRegression

#dt1 = pd.read_csv(os.path.join(savedir, 'NC_avgload_month.csv'), index_col=None)
# Loop over year-month to fit polynomial basis
# Version 1: standardize y in advance. Thus betas are effect times SD
dt0 = dt1.dropna(subset=['y_avg'])
a = -10.0; b = 40.0 # Common support for all months
mu_y = dt0['y_avg'].mean()
s_y = dt0['y_avg'].std()
dt0.loc[:, 'y_std'] = (dt0['y_avg']-mu_y)/s_y # standardize y to improve numerical stability

betas_list = []
months = dt0['yyyymm'].unique()
T = len(months)
for t in range(T):
    mon = months[t]
    dt_m = dt0.loc[dt0['yyyymm']==mon, ['yyyymm', 'y_std', 'temperature']]
    # Legendre polynomials on [-1, 1]
    r = dt_m['temperature'].values
    y_std = dt_m['y_std'].values
    r_nml = 2*(r-a)/(b-a)-1

    # Design matrix
    X_poly = legvander(r_nml, deg=99) # degrees from 0 to 99, m=100
    # Ridge regression to penalize large coefficients
    ridge = Ridge(alpha=1.0, fit_intercept=False) # debatable (3)
    model = ridge.fit(X_poly, y_std)
    betas = model.coef_ # Estimated coefficients
    betas_list.append(betas)

# Save the beta coefficients
betas_array = np.array(betas_list) # T x m
np.savetxt(os.path.join(savedir, 'list_ridge1betas100_month_ystd.csv'), betas_array, delimiter=',')
# End of script.