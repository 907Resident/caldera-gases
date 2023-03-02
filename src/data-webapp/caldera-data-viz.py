# %% import neccessary packages

# streamlit // web application package
import streamlit as st

# essential filepaths
from pathlib import Path
import pandas as pd
import numpy as np
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

from calderas import data_prep
from calderas import utils

importlib.reload(data_prep)
importlib.reload(utils)

# %% define page configuration

st.set_page_config(
    page_title="Caldera Gases",
    page_icon="🌎",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        'Report a bug': "https://github.com/907Resident/caldera-gases/issues",
        'About': "This is a web app for interacting with Caldera Gases."
    }
)

# %% Front matter

st.markdown(
    """
    # Caldera Gases

    A web application for viewing caldera gases data.

    ## Description

    To assist with the consumption of the work on caldera gases presented in this project, the authors have put together an interactive web application. This web application contains much of the data presented in the text related to this project and much more.

    The user of this application will gain the most value by reading through any prompts first before manipulating any of the figures. Additionally, there may be slight differences in the plots in the web app and publications due to randomness and/or slightly different data preparation steps. 
    
    """
)

# %% import data into workspace

st.markdown(
    """
    ## Preview the data

    View all available measurements collected by the researchers here in tabular form.
    """
)

# Boolean to resize the dataframe, stored as a session state variable
st.checkbox("Use container width", value=False, key="use_container_width")

main_df = utils.load_main_data()

# Display the dataframe and allow the user to stretch the dataframe
# across the full width of the container, based on the checkbox value
st.dataframe(main_df, use_container_width=st.session_state.use_container_width)
