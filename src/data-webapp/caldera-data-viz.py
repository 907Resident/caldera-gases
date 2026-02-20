# %% import necessary packages

import streamlit as st

# essential filepaths
from pathlib import Path
import os

# %% establish paths as variables

project_dir = Path(__file__).parents[2]
src_dir = Path(__file__).parents[1]
data_dir = os.path.join(project_dir, "pertinent-data")
