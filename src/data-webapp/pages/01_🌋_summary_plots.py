# %% import neccessary packages

# streamlit // web application package
import streamlit as st

# plotly // interactive visual package
from plotly.subplots import make_subplots

import plotly.graph_objects as go
import plotly.express as px

# urllib // package to interact with urls and websites
from urllib.request import urlopen

# pandas-geojson // converts traditional tabular data to geojson files
import pandas_geojson

# essential filepaths
from pathlib import Path

import pandas as pd
import numpy as np
import json
import os

# %% establish paths as variables

src_dir = Path(__file__).parents[1]
project_dir = Path(__file__).parents[2]
data_dir = os.path.join(project_dir, "pertinent-data")
webapp_dir = os.path.join(src_dir, "data-webapp")

# %% 

import importlib
import sys

sys.path.append(os.path.join(src_dir, "custom-package"))

from calderas import utils

importlib.reload(utils)

# %% 

main_df = utils.load_main_data()

main_df.dropna(inplace=True)

main_df = (
    main_df
    .rename(
        columns={
            "Latitude":"latitude",
            "Longitude":"longitude"
        }
    )
    .astype(
        {
            "latitude":float,
            "longitude":float
        }
    )
)

st.dataframe(main_df)

# %% 

geo_json = pandas_geojson.to_geojson(
    df=main_df,
    lat="latitude",
    lon="longitude",
    properties=[
        "Site_Name", "Soil_Classification", "CH4_Flux", "CO2_Flux"
    ]
)

# %% 

st.markdown(
    """
    ## Visualize the data in space

    """
)

st.map(
    main_df,
    zoom=5,
    use_container_width=True
)

# %% 

with urlopen('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_500k.json') as response:
    states = json.load(response)

# fig = px.choropleth(
#     main_df,
#     geojson=states,
#     locations="CH4_Flux",
#     color="Soil_Classification",
#     scope="usa"
# )

# fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
# st.plotly_chart(fig)

soil_type_color = {'Acid-Sulphate': '#cc33ff', 'Neutral-Chloride': '#00ccff', 'Unclassified': '#cc0000', 'Travertine': '#ffff00'}

main_df["plot_color_sclf"] = main_df.Soil_Classification.map(soil_type_color)

fig = go.Figure(
    go.Choroplethmapbox(
        geojson=states,
        locations=main_df.Site_Name,
        featureidkey="properties.name",
        z=main_df.CH4_Flux,
        colorscale="sunsetdark"
    )
)

fig.add_scattermapbox(
    lat=list(main_df.latitude),
    lon=list(main_df.longitude),
    mode = "markers+text",
    text=list(main_df.Site_Name),
    marker_size=10,
    marker=go.scattermapbox.Marker(
        color=list(main_df.plot_color_sclf)
    )
)

fig.update_layout(
    mapbox_style="carto-positron",
    mapbox_zoom=4.6,
    mapbox_center={"lat": 39.7392, "lon": -104.9903},
    width=800,
    height=600,
)

fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
st.plotly_chart(fig)