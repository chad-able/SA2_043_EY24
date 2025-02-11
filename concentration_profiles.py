import pandas as pd
import re
import init_functions
import geopandas as gpd

newts_filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
well_filename = "E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb"
shale_filename = "E:/codes/RC24/PWMapping/SedimentaryBasins_US_EIA/Lower_48_Sedimentary_Basins.shp"

gpd_args = {
    'layer': 1,
}

# initial data read
huc_data = init_functions.read_well_data(well_filename, **gpd_args)

# addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
huc_data = init_functions.centroids(huc_data, 0)

huc_data = init_functions.flow_rate(huc_data)

huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)

newts_data = init_functions.read_conc_data(newts_filename)

# concentration_profiles = init_functions.concentration_profiles_no_flow_rate(newts_data, shale_filename)

concentration_profiles_flow = init_functions.concentration_profiles_flow_rate(newts_data, shale_filename, huc_data)