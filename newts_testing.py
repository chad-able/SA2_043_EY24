import pandas as pd
import re
import init_functions
import geopandas as gpd

filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
df = init_functions.read_conc_data(filename)


basin_medians = df.groupby('BASIN')['Li'].median().reset_index()
field_medians = df.groupby('FIELD')['Li'].median().reset_index()
formation_medians = df.groupby('FORMATION')['Li'].median().reset_index()

basin_medians.columns = ['name', 'Li_concentration']
field_medians.columns = ['name', 'Li_concentration']
formation_medians.columns = ['name', 'Li_concentration']

basin_medians['type'] = 'basin'
formation_medians['type'] = 'formation'
field_medians['type'] = 'field'

all_medians = pd.concat([basin_medians, field_medians, formation_medians], ignore_index=True)

all_medians.columns = ['name', 'Li_concentration', 'type']

all_medians.to_csv('Li_medians.csv', index=False)

plays = ['Permian']

df = init_functions.filter_conc_data(df, plays)

df.to_csv('test.csv')

usgs_shapefile = "E:/codes/RC24/PWMapping/USGS_NPWGDv3_shape/USGS_NPWGDv3_shape.shp"

well_filename = "E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb"

gpd_args = {
    'layer': 0,
    'rows': 5000
}

shale_plays = gpd.read_file(usgs_shapefile, **gpd_args)

print(shale_plays.columns)

df.to_csv('test.csv')
# gpd_args = {
#     'layer': 1,
#     'rows': 500
# }
#
#
# # initial data read
# huc_data = init_functions.read_well_data(well_filename, **gpd_args)
#
# #addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
# huc_data = init_functions.centroids(huc_data, 0)
#
#
#
#
#
#
# huc_data = init_functions.flow_rate(huc_data)
#
# huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)
#
#
#
#
# huc_data = init_functions.shale_plays(huc_data, well_filename)
#
# print(huc_data.head())