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
savedir = os.path.abspath(os.path.join(thisdir, '..', 'output', 'month'))
os.makedirs(savedir, exist_ok=True)
#%%
dt1 = pd.read_csv(os.path.join(savedir, 'NC_avgload_month.csv'), index_col=None) # change here for different version
dt1['datetime_UTC'] = pd.to_datetime(dt1['datetime_UTC'])
dt1['Year'] = dt1['Year'].astype(int)
dt1['Month'] = dt1['Month'].astype(int)
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

# Plot a sorted line for each month
for i, mon in enumerate(months):
    dt_m = dt1[dt1['yyyymm'] == mon].copy()
    # Sort by temperature
    dt_m = dt_m.sort_values('temperature')
    # Get month index for color
    month_idx = int(dt_m['Month'].iloc[0]) - 1  # 0-based
    color_rgb = cmap(month_idx)[:3]
    color_str = f'rgb({int(color_rgb[0]*255)}, {int(color_rgb[1]*255)}, {int(color_rgb[2]*255)})'
    fig.add_trace(go.Scatter3d(
        x=dt_m['temperature'],
        y=[mon]*len(dt_m),
        z=dt_m['y_avg'],
        mode='lines',
        line=dict(color=color_str, width=2),
        name=calendar.month_abbr[month_idx+1] if i < 12 else str(mon),
        showlegend=False,
        hoverinfo='skip'
    ))

# Manually add 12 dummy traces for legend
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
years = sorted(set(dt1['Year']))
fig.update_layout(
    scene=dict(
        xaxis_title='Temperature (C)',
        yaxis=dict(
            title='Year',
            tickvals=[dt1[dt1['Year'] == y]['yyyymm'].iloc[0] for y in years],
            ticktext=[str(y) for y in years]
        ),        
        zaxis_title='Load (MW)'
    ),
    legend_title_text='Month',
    margin=dict(l=0, r=0, b=0, t=30),
    height=800,
    title='3D Plot: Temperature vs Smoothed Load by Month'
)

fig.write_html(os.path.join(savedir, 'interactive_smooth_plot_2.html')) # change the name for different version
# End of script.