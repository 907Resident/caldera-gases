# %% import necessary libraries and packages


# essential packages
import pandas as pd
import os

# %% define function for basic preprocessing

def basic_preprocessing(df):
    """
    The basic_preprocessing function takes a dataframe as an argument and performs the following:
        1. Set dictionary for correct data types
        2. Set data types accordingly
        3. Drop unnecessary columns
    
    :param df: Pass the data frame to be preprocessed
    :return: The data frame with the following columns:
    :doc-author: Trelent & Moyo Ajayi
    """
    
    # Set dictionary for correct data types
    data_types = {
        "Group":"category",
        "Location":"object",
        "Soil_Classification":"category",
        "Latitude":"float64","Longitude":"float64",
        "Date_of_Measurement":"datetime64",
        "Start_Time_of_Chamber_Enclosure":"datetime64",
        "End_Time_of_Chamber_Enclosure":"datetime64",
        "Duration_of_Total_Chamber_Enclosure":"float64",
        "Start_Time_Flux":"datetime64",
        "End_Time_Flux":"datetime64",              
        "CH4_Flux":"float64",
        "LowerBound_CH4_Flux":"float64",
        "UpperBound_CH4_Flux":"float64",
        "KP_d13CH4_source":"float64",
        "KP_LowerBound_d13CH4_source":"float64",
        "KP_UpperBound_d13CH4_source":"float64",
        "d13CH4_source":"float64", 
        "d13CH4_LowerBound_source":"float64",
        "d13CH4_UpperBound_source":"float64",
        "CO2_Flux":"float64",
        "LowerBound_CO2_Flux":"float64",
        "UpperBound_CO2_Flux":"float64",
        "d13CO2_source":"float64",
        "KP_d13CO2_source":"float64",
        "KP_LowerBound_d13CO2_source":"float64",
        "KP_UpperBound_d13CO2_source":"float64",
        "d13CO2_source":"float64",
        "d13CO2_LowerBound_source":"float64",
        "d13CO2_UpperBound_source":"float64",
        "Horita_Geothermometer_Temperature_at_Formation":"float64",
        "Ambient_Temperature":"float64",
        "Barometric_Pressure":"float64",
        "Soil_Tempeature_at_Surface":"float64"
    }
    # Set data types accordingly
    df = df.astype(data_types)
    # Drop unnecessary columns
    df.drop(df.columns[30:33], axis=1, inplace=True)
    df.drop(df.columns[-1], axis=1, inplace=True)

    # Return data frame
    return df