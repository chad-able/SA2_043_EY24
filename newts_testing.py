import pandas as pd
import re
import init_functions
import geopandas as gpd

filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
df = init_functions.read_conc_data(filename)

newts_usgs = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df['LONGITUDE'], df['LATITUDE']),crs="EPSG:4269")
# #
# # basin_medians = df.groupby('BASIN')['Li'].median().reset_index()
# # field_medians = df.groupby('FIELD')['Li'].median().reset_index()
# # formation_medians = df.groupby('FORMATION')['Li'].median().reset_index()
# #
# # basin_medians.columns = ['name', 'Li_concentration']
# # field_medians.columns = ['name', 'Li_concentration']
# # formation_medians.columns = ['name', 'Li_concentration']
# #
# # basin_medians['type'] = 'basin'
# # formation_medians['type'] = 'formation'
# # field_medians['type'] = 'field'
# #
# # all_medians = pd.concat([basin_medians, field_medians, formation_medians], ignore_index=True)
# #
# # all_medians.columns = ['name', 'Li_concentration', 'type']
# #
# # all_medians.to_csv('Li_medians.csv', index=False)
# #
# # plays = ['Permian']
# #
# # df = init_functions.filter_conc_data(df, plays)
# #
# # df.to_csv('test.csv')
#
# # usgs_shapefile = "E:/codes/RC24/PWMapping/USGS_NPWGDv3_shape/USGS_NPWGDv3_shape.shp"
# #
# # usgs_data = init_functions.read_conc_shapefile(usgs_shapefile)
# #
# # print(usgs_data.crs)
#
well_filename = "E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb"

gpd_args = {
    'layer': 1,
}

# # shale_plays = gpd.read_file(usgs_shapefile, **gpd_args)
# #
# # print(shale_plays.columns)
# #
# # df.to_csv('test.csv')
# #
# #
# # initial data read
huc_data = init_functions.read_well_data(well_filename, **gpd_args)

#addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
huc_data = init_functions.centroids(huc_data, 0)

huc_data = init_functions.flow_rate(huc_data)

huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)
# #
# # usgs_data.to_crs(epsg=5070)
# # usgs_data['centroid'] = usgs_data.geometry.centroid
# #
# # usgs_centroids = usgs_data.copy()
# # usgs_centroids.set_geometry('centroid', inplace=True)
# # usgs_data = usgs_data.to_crs(huc_data.crs)
# newts_usgs = newts_usgs.to_crs(huc_data.crs)
# #
# # print(usgs_data[usgs_data['Li'] == -9999])
# # print(usgs_data[~usgs_data.is_valid])
# #

conc_array = ['TDS_combined', 'Ag', 'Al', 'As', 'Au', 'B', 'BO3', 'Ba', 'Be', 'Bi', 'Br', 'CO3', 'HCO3', 'Ca', 'Cd', 'Cl', 'Co',
              'Cr', 'Cs', 'Cu', 'F', 'FeTot', 'FeIII', 'FeII', 'FeS', 'FeAl', 'FeAl2O3', 'Hg', 'I', 'K', 'KNa', 'Li', 'Mg',
              'Mn', 'Mo', 'N', 'NO2', 'NO3', 'NO3NO2', 'NH4', 'TKN', 'Na', 'Ni', 'OH', 'P', 'PO4', 'Pb', 'Rh', 'Rb', 'S', 'SO3',
              'SO4', 'HS', 'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ti', 'Tl', 'U', 'V', 'W', 'Zn', 'ALKHCO3', 'ACIDITY', 'DIC', 'DOC',
              'TOC', 'CN', 'BOD', 'COD', 'BENZENE', 'TOLUENE', 'ETHYLBENZ', 'XYLENE', 'ACETATE', 'BUTYRATE', 'FORMATE', 'LACTATE', 'PHENOLS',
              'PERC', 'PROPIONATE', 'PYRUVATE', 'VALERATE', 'ORGACIDS', 'Ar', 'CH4', 'C2H6', 'CO2', 'H2', 'H2S', 'He', 'N2', 'NH3',
              'O2', 'ALPHA', 'BETA', 'dD', 'H3', 'd7Li', 'd11B', 'd13C', 'C14', 'd18O', 'd34S', 'd37Cl', 'K40', 'd81Br',
              'Sr87Sr86', 'I129', 'Rn222', 'Ra226', 'Ra228']

joined = gpd.sjoin(newts_usgs, huc_data, how="inner", predicate="within")
median_concentration = joined.groupby("huc8")[conc_array].median().reset_index()
# #
huc_data = huc_data.merge(median_concentration, on="huc8", how="left")

huc_data[['huc8', 'Li']].to_csv('huc8_Li.csv', index=False)
#
shale_filename = "E:/codes/RC24/PWMapping/SedimentaryBasins_US_EIA/Lower_48_Sedimentary_Basins.shp"






# shale_plays = gpd.read_file(shale_filename)
#
# for i in range(len(shale_plays.columns)):
#     print(shale_plays.columns[i])
#
#
# #
# #
# #
# #
huc_data = init_functions.shale_plays(huc_data, shale_filename)
#
plays = ['Permian']
# fix conc data calculations
# newts_data = init_functions.filter_conc_data(newts_data, plays)
# print(newts_data['Li'].median())
pattern = '|'.join(re.escape(play) for play in plays)

huc_data = huc_data[huc_data['Shale_play'].str.contains(pattern, case=False, na=False)]

for i in range(len(conc_array)):

    valid = huc_data[huc_data[conc_array[i]].notnull()]

    if len(valid) > 0:

        weighted_avg = (valid[conc_array[i]] * valid['2022_flow_gpm']).sum() / valid['2022_flow_gpm'].sum()

        huc_data[conc_array[i]] = huc_data[conc_array[i]].fillna(weighted_avg)

huc_data[['huc8', 'Li', 'Shale_play']].to_csv('huc8_Li.csv', index=False)
#
# print(huc_data.head())