# %% import necessary libraries and packages

from pathlib import Path
import streamlit as st

# essential packages
import pandas as pd
import os

# %% 

project_dir = Path(__file__).parents[3]
data_dir = os.path.join(project_dir, "pertinent-data")

# %% define function to load data

@st.cache_data
def load_main_data():
    return pd.read_csv(
        os.path.join(data_dir, "summary-data", "all_data_cleaned.csv")
    )