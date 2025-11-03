#%%
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.colors import Normalize
import plotly.graph_objects as go
import calendar  # to get month names
###########################
### Interactive 3D plot ### 
# Color gradient by month #
#%%
thisdir = os.path.dirname(__file__)
savedir = os.path.abspath(os.path.join(thisdir, '..', 'output'))
os.makedirs(savedir, exist_ok=True)
#%%
url = "https://raw.githubusercontent.com/BrianIZKom1911/electricload_CTRF/master/data_clean/NC_main.csv"
data = pd.read_csv(url, index_col=None)
data['datetime_UTC'] = pd.to_datetime(data['datetime_UTC'])
data['Year'] = data['Year'].astype(int)
data['Month'] = data['Month'].astype(int)
data['yyyymm'] = data['Year'].astype(str) + '-' + data['Month'].astype(str).str.zfill(2)
#data['logy'] = np.log(data['load'])
dt1 = data[(data['temperature'] >= -10) & (data['temperature'] <= 40)]
#%%
# Set up
cmap = mpl.colormaps['turbo'].resampled(12)
norm = Normalize(vmin=0, vmax=12)
colors = [cmap(i/12) for i in range(12)]
month_order = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1] # Custom order for months
custom_colors = [colors[i] for i in month_order]
custom_colorscale = [
    [i/11, f'rgb({int(c[0]*255)}, {int(c[1]*255)}, {int(c[2]*255)})'] 
    for i, c in enumerate(custom_colors)
]

fig = go.Figure()
months = dt1['yyyymm'].unique()

fig.add_trace(go.Scatter3d(
    x=dt1['temperature'],
    y=dt1['yyyymm'].astype('category').cat.codes, # convert dates to numeric codes
    z=dt1['load'],
    mode='markers',
    marker=dict(
        size=1,
        color=dt1['Month'], # Color by month
        colorscale=custom_colorscale,
        cmin=1, # January
        cmax=12 # December
        #opacity=0.7,
        #colorbar=dict(title='Month')
    ),
    hovertext=dt1['yyyymm'],  # show date on hover
    hoverinfo='x+z+text'
))

# Manually add 12 dummy traces for legend # to be consistent with line plot
for i in range(12):
    color_rgb = cmap(i)[:3]
    color_str = f'rgb({int(color_rgb[0]*255)}, {int(color_rgb[1]*255)}, {int(color_rgb[2]*255)})'
    fig.add_trace(go.Scatter3d(
        x=[None], y=[None], z=[None],
        mode='lines',
        line=dict(color=color_str, width=4),
        name=calendar.month_abbr[i+1], # Jan, Feb, ... Dec
        showlegend=True
    ))

# Customize axes
fig.update_layout(
    scene=dict(
        xaxis_title='Temperature (C)',
        yaxis=dict(
            title='Year-Month',
            tickvals=[0, 60, 120, 180, 240],
            ticktext=['2002', '2007', '2012', '2017', '2022']
        ),        
        zaxis_title='Load (MW)'
    ),
    legend_title_text='Month',
    margin=dict(l=0, r=0, b=0, t=30),
    height=800,
    title='3D Scatterplot: Temperature vs Load by Month',
)

fig.write_html(os.path.join(savedir, 'interactive_scatterplot_1.html'))
# End of script.