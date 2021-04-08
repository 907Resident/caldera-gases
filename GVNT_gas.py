
#%% Establish basic libraries
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

#%%  Import accumulation chamber data
df_AC = pd.read_csv("E:\moyoa\Documents\OneDrive - Vanderbilt\Python\Yellowstone-Python\GVNT\GVNT_Lewicki\Lewicki_Soil_CO2_flux+temperature_GVNT.csv",
parse_dates=True, verbose=True)

# Delete second row
df_AC.drop([0], axis=0, inplace=True)
# Reset Index
df_AC.reset_index(inplace=True, drop=True)

# Change datatypes
## Date
df_AC["Date"] = df_AC.Date.astype("datetime64")
## Numeric Data
df_AC.iloc[:,1:] = df_AC.iloc[:,1:].apply(pd.to_numeric)

# Check Data Types
df_AC.dtypes

#%% Convert UTM measurements to Lat/Lon
# Using the `utm` package, convert the data from UTM to Lat/Lon (decimal). This will ensure they are easier 
# to plot when mapping. Function below grabbed from [Stack](https://stackoverflow.com/questions/49890492/convert-utm-to-lat-long-in-csv-using-pandas). 
#%% utm2latlon - function
import utm
# Convert UTM to Lat/Lon
def utm2latlon(df, zone_num, zone_letter):
    '''
    df = dataframe conatining the UTM measurements with headers for "Easting" and "Norting"
    zone_num = UTM Zone Number
    zone_letter = UTM Zone Letter

    Returns df with columns inserted into df as "lat" and "lon"
    '''
    lat, lon = utm.to_latlon(df["Easting"], df["Northing"], zone_number=zone_num, zone_letter=zone_letter)
    lat = lat.astype("float")
    lon = lon.astype("float")
    df["lat"] = lat
    df["lon"] = lon
    
    # Return dataframe
    return df

#%%  Apply the UTM Conversion to the data
df_AC = utm2latlon(df_AC, 12, "T")
df_AC.head()

#%% Quick Charts
# Some charts to view the fluxes (g m-2 day-1) and temperature (°C)

## Filter out -9999.0 as NaN
df_AC = df_AC.replace(-9999.0, np.NaN)

# Set Theme (e.g., visual style, color pallette, etc.)
sns.set_theme(context="paper", style="darkgrid", palette="rainbow",
              font_scale=1.5)
# Establish path to save figure
fig_path = "E:\moyoa\Documents\OneDrive - Vanderbilt\Python\Yellowstone-Python\GVNT\GVNT_Lewicki"
# KDE of Fluxes
sns.kdeplot(data=df_AC, x="Soil CO2 flux", palette="rainbow", fill=True)
plt.show()
plt.savefig(fig_path+"\Lewicki_GVNT_CO2_AC_fluxes_kde.png", dpi=300)

# KDE of Temperatures
sns.kdeplot(data=df_AC, x="Soil temperature at 30 cm depth", palette="rainbow", fill=True)
plt.show()
plt.savefig(fig_path+"\Lewicki_GVNT_SoilTemps_kde.png", dpi=300)

# Scatterplot of Fluxes and Temperature
sns.scatterplot(data=df_AC, x="Soil temperature at 30 cm depth", y="Soil CO2 flux")
plt.yscale("log")
plt.ylabel("CO2 Flux (g m-2 day-1)")
plt.xlabel("Soil Temp (°C)")
plt.show()
plt.savefig(fig_path+"\Lewicki_GVNT_CO2_flux_v_SoilTemp_scatt.png", dpi=300)

#%% Plot Data in Space
# Using the [`geopandas`](https://geopandas.org/) module, we can plot Lewicki's data on a simple map. 
# Our data will be  superimposed on this map as well for comparison.

# Note: There was a ton of trouble to get the `geopandas` package to load properly. See this dicussion 
# (https://github.com/Toblerity/Fiona/issues/944) for more details
#%% Import necessary libraries
import fiona
import geopandas as gpd
import contextily as ctx 
from shapely.geometry import Point, Polygon
#%% Establish a geopandas dataframe
# As described in this tutorial, one can convert a pandas dataframe into a geopandas dataframe by following
# the steps presented [here](https://medium.com/@ianforrest11/graphing-latitudes-and-longitudes-on-a-map-bf64d5fca391).
#%% 

# Import Caldera Boundaries
caldera_shp_path = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\YNP Maps\Caldera\calderas.shp"
caldera_boundary = gpd.read_file(caldera_shp_path)

# Designate coordinate system
crs = {"init":"epsg:4326"}

# Zip x and y coordinates into single feature
geometry = [Point(xy) for xy in zip(df_AC["lon"], df_AC["lat"])]

# Create GeoPandas dataframe
geo_df = gpd.GeoDataFrame(df_AC, crs = crs, geometry = geometry)

# Drop UTM
geo_df.drop(columns=["Easting", "Northing"], inplace=True)

# Preview dataframe
geo_df.head()

#%% 
# Import Tick Label formatter
from matplotlib.ticker import FormatStrFormatter

# Plot
fig, ax = plt.subplots(figsize=(10,10))

# Add .shp file to map
caldera_boundary.plot(ax=ax, alpha=0.6)
ctx.add_basemap(ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.3, zoom=5)

# Add data
geo_df.plot(column="Soil CO2 flux", ax=ax, alpha=0.8, legend=True, markersize=10, cmap="nipy_spectral_r")
ax.xaxis.set_major_formatter(FormatStrFormatter('%3.4f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%3.4f'))
# Set latitiude and longitude boundaries for map display
plt.xlim(geo_df.lon.min(), geo_df.lon.max())
plt.ylim(geo_df.lat.min(), geo_df.lat.max())
plt.show()

# %%
