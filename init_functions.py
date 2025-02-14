import pandas as pd
import geopandas as gpd
import random
import pyomo.environ as pyo
import modeling_functions
import re
import print_functions


# TODO
# add functionality to separate out single play, list of plays
# drop wells that do not meet certain flow rate requirements
# Add in NEWTS concentration data, match to shale play
# Fix shale play relations

def read_well_data(filename, **gpd_args):
    gpdf = gpd.read_file(filename, **gpd_args)
    return gpdf

def read_conc_data(filename):
    df = pd.read_csv(filename, encoding='ISO-8859-1', low_memory = False)
    newts_usgs = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['LONGITUDE'], df['LATITUDE']), crs="EPSG:4269")
    return newts_usgs

def read_conc_shapefile(filename, **gpd_args):
    df = gpd.read_file(filename, **gpd_args)
    return df

def filter_conc_data(df, plays):
    pattern = '|'.join(re.escape(play) for play in plays)

    df = df[df[['BASIN', 'FIELD', 'FORMATION']].apply(
        lambda row: row.str.contains(pattern, case=False, na=False).any(), axis=1
    )]

    return df


def combine_conc_data(huc_df, conc_df):
    conc_array = ['TDS_combined', 'Ag', 'Al', 'As', 'Au', 'B', 'BO3', 'Ba', 'Be', 'Bi', 'Br', 'CO3', 'HCO3', 'Ca', 'Cd',
                  'Cl', 'Co',
                  'Cr', 'Cs', 'Cu', 'F', 'FeTot', 'FeIII', 'FeII', 'FeS', 'FeAl', 'FeAl2O3', 'Hg', 'I', 'K', 'KNa',
                  'Li', 'Mg',
                  'Mn', 'Mo', 'N', 'NO2', 'NO3', 'NO3NO2', 'NH4', 'TKN', 'Na', 'Ni', 'OH', 'P', 'PO4', 'Pb', 'Rh', 'Rb',
                  'S', 'SO3',
                  'SO4', 'HS', 'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ti', 'Tl', 'U', 'V', 'W', 'Zn', 'ALKHCO3',
                  'ACIDITY', 'DIC', 'DOC',
                  'TOC', 'CN', 'BOD', 'COD', 'BENZENE', 'TOLUENE', 'ETHYLBENZ', 'XYLENE', 'ACETATE', 'BUTYRATE',
                  'FORMATE', 'LACTATE', 'PHENOLS',
                  'PERC', 'PROPIONATE', 'PYRUVATE', 'VALERATE', 'ORGACIDS', 'Ar', 'CH4', 'C2H6', 'CO2', 'H2', 'H2S',
                  'He', 'N2', 'NH3',
                  'O2', 'ALPHA', 'BETA', 'dD', 'H3', 'd7Li', 'd11B', 'd13C', 'C14', 'd18O', 'd34S', 'd37Cl', 'K40',
                  'd81Br',
                  'Sr87Sr86', 'I129', 'Rn222', 'Ra226', 'Ra228']

    conc_df = conc_df.to_crs(huc_df.crs)
    joined = gpd.sjoin(conc_df, huc_df, how="inner", predicate="within")
    median_concentration = joined.groupby("huc8")[conc_array].median().reset_index()
    # #
    huc_df = huc_df.merge(median_concentration, on="huc8", how="left")
    return huc_df


def centroids(gpdf, proj_flag):

    # set up for multiple projections - distance equations are currently only set up for EPSG 4326

    if proj_flag == 0:
        crs = gpdf["geometry"].to_crs(epsg=4326)
    elif proj_flag == 1:
        crs = huc_data["geometry"].to_crs('+proj=cea')
    centroids = crs.centroid
    # convert back
    # centroids = centroids.to_crs(huc8_data.crs)
    gpdf["geometry_centroids"] = centroids
    gpdf['lat'] = centroids.y
    gpdf['lon'] = centroids.x
    return gpdf

def flow_rate(gpdf, remove_zeros = True, remove_nans = True):

    # built specifically for the well summary file
    gpdf["2022_flow_gpm"] = gpdf["F2022_WaterProd_BBL_sum_1"].apply(flow_rate_calc)

    if remove_zeros:
        gpdf = gpdf[(gpdf["2022_flow_gpm"] != 0)]

    if remove_nans:
        gpdf = gpdf[(gpdf["2022_flow_gpm"].notna())]

    return gpdf

def flow_rate_calc(x):
    #converts bbl/year to gal/min
    x_per_min = x / 365 / 24 / 60
    x_gal_per_min = x_per_min / 42
    return x_gal_per_min

def shale_plays(gpdf, filename):
    gpd_args = {
        'layer': 0,
    }
    shale_plays = gpd.read_file(filename, **gpd_args)
    shale_plays.rename(columns={'Name': 'Shale_play'}, inplace=True)
    #save original CRS
    originalcrs = gpdf.crs

    #change to shale CRS
    water_locations_gdf = gpdf.to_crs(shale_plays.crs)


    water_locations_gdf["geometry"] = water_locations_gdf.buffer(0)
    shale_plays["geometry"] = shale_plays.buffer(0)
    validity_check(water_locations_gdf)
    validity_check(shale_plays)
    #join dataframes


    water_in_shale = gpd.sjoin(water_locations_gdf, shale_plays, how="left", predicate="intersects")
    print(water_in_shale.index)
    print(water_locations_gdf.index)
    print(shale_plays['Shale_play'])

    water_in_shale = water_in_shale.reset_index()
    # point_summary = intersection_percentage(water_in_shale, shale_plays)
    # print(point_summary)
    # Group by water location index and aggregate shale play names
    water_in_shale_grouped = water_in_shale.groupby("index").agg({
        'Shale_play': lambda x: list(set(x.dropna()))  # Collect non-NaN matches as a list
    }).reset_index()

    # Merge back to the original dataset
    water_locations_with_shales = water_locations_gdf.merge(water_in_shale_grouped, left_index=True, right_on='index',
                                                            how="left")

    # Optional: Convert list to string for better readability
    water_locations_with_shales['Shale_play'] = water_locations_with_shales['Shale_play'].apply(
        lambda x: ', '.join(x) if isinstance(x, list) else None
    )

    print(water_locations_with_shales.head())
    #restore original CRS
    restored_gpdf = water_locations_with_shales.to_crs(originalcrs)

    # water_locations_with_shales.to_csv("water_locations.csv", index=False)

    return restored_gpdf

def concentration_profiles_no_flow_rate(newts_df, filename):

    conc_array = ['TDS_combined', 'Ag', 'Al', 'As', 'Au', 'B', 'BO3', 'Ba', 'Be', 'Bi', 'Br', 'CO3', 'HCO3', 'Ca', 'Cd',
                  'Cl', 'Co',
                  'Cr', 'Cs', 'Cu', 'F', 'FeTot', 'FeIII', 'FeII', 'FeS', 'FeAl', 'FeAl2O3', 'Hg', 'I', 'K', 'KNa',
                  'Li', 'Mg',
                  'Mn', 'Mo', 'N', 'NO2', 'NO3', 'NO3NO2', 'NH4', 'TKN', 'Na', 'Ni', 'OH', 'P', 'PO4', 'Pb', 'Rh', 'Rb',
                  'S', 'SO3',
                  'SO4', 'HS', 'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ti', 'Tl', 'U', 'V', 'W', 'Zn', 'ALKHCO3',
                  'ACIDITY', 'DIC', 'DOC',
                  'TOC', 'CN', 'BOD', 'COD', 'BENZENE', 'TOLUENE', 'ETHYLBENZ', 'XYLENE', 'ACETATE', 'BUTYRATE',
                  'FORMATE', 'LACTATE', 'PHENOLS',
                  'PERC', 'PROPIONATE', 'PYRUVATE', 'VALERATE', 'ORGACIDS', 'Ar', 'CH4', 'C2H6', 'CO2', 'H2', 'H2S',
                  'He', 'N2', 'NH3',
                  'O2', 'ALPHA', 'BETA', 'dD', 'H3', 'd7Li', 'd11B', 'd13C', 'C14', 'd18O', 'd34S', 'd37Cl', 'K40',
                  'd81Br',
                  'Sr87Sr86', 'I129', 'Rn222', 'Ra226', 'Ra228']

    gpd_args = {
        'layer': 0,
    }
    shale_df = gpd.read_file(filename, **gpd_args)
    shale_df.rename(columns={'Name': 'Shale_play'}, inplace=True)
    newts_df = newts_df.to_crs(shale_df.crs)
    joined = gpd.sjoin(newts_df, shale_df, how="inner", predicate="within")
    median_concentration = joined.groupby("Shale_play")[conc_array].median().reset_index()
    # #
    shale_df = shale_df.merge(median_concentration, on="Shale_play", how="left")

    shale_array = ['Shale_play'] + conc_array
    shale_df[shale_array].to_csv('concentration_profiles_no_flow.csv', index=False)

    return shale_df

def concentration_profiles_flow_rate(newts_df, shale_filename, well_df):

    conc_array = ['TDS_combined', 'Ag', 'Al', 'As', 'Au', 'B', 'BO3', 'Ba', 'Be', 'Bi', 'Br', 'CO3', 'HCO3', 'Ca', 'Cd',
                  'Cl', 'Co',
                  'Cr', 'Cs', 'Cu', 'F', 'FeTot', 'FeIII', 'FeII', 'FeS', 'FeAl', 'FeAl2O3', 'Hg', 'I', 'K', 'KNa',
                  'Li', 'Mg',
                  'Mn', 'Mo', 'N', 'NO2', 'NO3', 'NO3NO2', 'NH4', 'TKN', 'Na', 'Ni', 'OH', 'P', 'PO4', 'Pb', 'Rh', 'Rb',
                  'S', 'SO3',
                  'SO4', 'HS', 'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ti', 'Tl', 'U', 'V', 'W', 'Zn', 'ALKHCO3',
                  'ACIDITY', 'DIC', 'DOC',
                  'TOC', 'CN', 'BOD', 'COD', 'BENZENE', 'TOLUENE', 'ETHYLBENZ', 'XYLENE', 'ACETATE', 'BUTYRATE',
                  'FORMATE', 'LACTATE', 'PHENOLS',
                  'PERC', 'PROPIONATE', 'PYRUVATE', 'VALERATE', 'ORGACIDS', 'Ar', 'CH4', 'C2H6', 'CO2', 'H2', 'H2S',
                  'He', 'N2', 'NH3',
                  'O2', 'ALPHA', 'BETA', 'dD', 'H3', 'd7Li', 'd11B', 'd13C', 'C14', 'd18O', 'd34S', 'd37Cl', 'K40',
                  'd81Br',
                  'Sr87Sr86', 'I129', 'Rn222', 'Ra226', 'Ra228']

    gpd_args = {
        'layer': 0,
    }
    shale_df = gpd.read_file(shale_filename, **gpd_args)
    shale_df.rename(columns={'Name': 'Shale_play'}, inplace=True)
    newts_df = newts_df.to_crs(shale_df.crs)
    joined = gpd.sjoin(newts_df, shale_df, how="inner", predicate="within")
    median_concentration = joined.groupby("Shale_play")[conc_array].median().reset_index()
    # #

    well_df = well_df.to_crs(newts_df.crs)
    well_intersection = gpd.overlay(well_df, shale_df, how="intersection")
    well_intersection["intersect_area"] = well_intersection.geometry.area
    well_df["orig_area"] = well_df.geometry.area
    well_intersection = well_intersection.merge(well_df[["huc8", "orig_area", "2022_flow_gpm"]], on="huc8", how="left")
    well_intersection.rename(columns={'2022_flow_gpm_x': '2022_flow_gpm', '2022_flow_gpm_y': '2022_flow_gpm'}, inplace=True)
    well_intersection = well_intersection.loc[:, ~well_intersection.columns.duplicated()]
    well_intersection["weight"] = well_intersection["intersect_area"] / well_intersection["orig_area"]
    for i in range(len(well_intersection.columns)):
        print(well_intersection.columns[i])
    well_intersection["weighted_flow"] = well_intersection["2022_flow_gpm"] * well_intersection["weight"]
    flow_stats = well_intersection.groupby("Shale_play").agg(
        total_weighted_flow=("weighted_flow", "sum"),
        median_flow=("2022_flow_gpm", "median")
    ).reset_index()

    basin_stats = median_concentration.merge(flow_stats, on="Shale_play", how="left")

    shale_array = ['Shale_play', 'total_weighted_flow', 'median_flow'] + conc_array
    basin_stats[shale_array].to_csv('concentration_profiles_flow.csv', index=False)
    print_functions.mapping_flow_rates(shale_df, well_df, basin_stats)
    return basin_stats


def flow_rate_dump(gpdf):
    return gpdf

def filter_by_shale_plays(df, plays, replace_values = True):
    pattern = '|'.join(re.escape(play) for play in plays)

    df = df[df['Shale_play'].str.contains(pattern, case=False, na=False)]

    conc_array = ['TDS_combined', 'Ag', 'Al', 'As', 'Au', 'B', 'BO3', 'Ba', 'Be', 'Bi', 'Br', 'CO3', 'HCO3', 'Ca', 'Cd',
                  'Cl', 'Co',
                  'Cr', 'Cs', 'Cu', 'F', 'FeTot', 'FeIII', 'FeII', 'FeS', 'FeAl', 'FeAl2O3', 'Hg', 'I', 'K', 'KNa',
                  'Li', 'Mg',
                  'Mn', 'Mo', 'N', 'NO2', 'NO3', 'NO3NO2', 'NH4', 'TKN', 'Na', 'Ni', 'OH', 'P', 'PO4', 'Pb', 'Rh', 'Rb',
                  'S', 'SO3',
                  'SO4', 'HS', 'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ti', 'Tl', 'U', 'V', 'W', 'Zn', 'ALKHCO3',
                  'ACIDITY', 'DIC', 'DOC',
                  'TOC', 'CN', 'BOD', 'COD', 'BENZENE', 'TOLUENE', 'ETHYLBENZ', 'XYLENE', 'ACETATE', 'BUTYRATE',
                  'FORMATE', 'LACTATE', 'PHENOLS',
                  'PERC', 'PROPIONATE', 'PYRUVATE', 'VALERATE', 'ORGACIDS', 'Ar', 'CH4', 'C2H6', 'CO2', 'H2', 'H2S',
                  'He', 'N2', 'NH3',
                  'O2', 'ALPHA', 'BETA', 'dD', 'H3', 'd7Li', 'd11B', 'd13C', 'C14', 'd18O', 'd34S', 'd37Cl', 'K40',
                  'd81Br',
                  'Sr87Sr86', 'I129', 'Rn222', 'Ra226', 'Ra228']

    if replace_values:
        for i in range(len(conc_array)):

            valid = df[df[conc_array[i]].notnull()]

            if len(valid) > 0:
                weighted_avg = (valid[conc_array[i]] * valid['2022_flow_gpm']).sum() / valid['2022_flow_gpm'].sum()

                df[conc_array[i]] = df[conc_array[i]].fillna(weighted_avg)

    return df

def validity_check(gdf):
    invalid_locations = gdf[~gdf.is_valid]
    print(f"Invalid water locations: {len(invalid_locations)}")

def intersection_percentage(gdf, shale_plays):

    intersection = gpd.overlay(gdf, shale_plays, how="intersection")
    intersection["intersection_area"] = intersection.geometry.area
    intersection.rename(columns={'Shale_play_1': 'Shale_play', 'Shale_play_2': 'Shale_play'}, inplace=True)
    intersection['Shale_play'] = intersection['Shale_play'].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
    )
    intersection['Shale_play'].fillna('Unknown', inplace=True)

    intersection = intersection.loc[:, ~intersection.columns.duplicated()]

    point_summary = intersection.groupby(['index', 'Shale_play']).agg({
        'intersection_area': 'sum'
    }).reset_index()

    buffer_areas = gdf.set_index('index').geometry.area
    buffer_areas = buffer_areas[~buffer_areas.index.duplicated()]
    point_summary['percentage'] = (
            point_summary['intersection_area'] / point_summary['index'].map(buffer_areas).fillna(0)) * 100

    point_summary.to_csv("Point_summary.csv", index=False)
    return point_summary


def model_init_any_location(lat_min, lat_max, lon_min, lon_max, num_facilities, num_sites, is_hybrid):
    model = pyo.ConcreteModel()

    def random_lat_init(model, i):
        return random.uniform(lat_min, lat_max)

    def random_lon_init(model, i):
        return random.uniform(lon_min, lon_max)

    def init_x_wrapper(model, i, j):
        return modeling_functions.initialize_x(model, i, j, num_facilities=num_facilities, num_sites=num_sites)

    def init_y_wrapper(model, j):
        return modeling_functions.initialize_y(model, j, num_sites = num_sites)
    def coverage_rule(model):
        return sum(model.x[i, j] for i in range(num_facilities) for j in range(num_sites)) >= num_sites

    # for centralized facilities (solids or liquids)
    model.x = pyo.Var(range(num_facilities), range(num_sites), domain=pyo.Binary, initialize=init_x_wrapper)
    if is_hybrid:
        # for deciding centralized or modular)
        model.y = pyo.Var(range(num_sites), domain=pyo.Binary, initialize=init_y_wrapper)
        # second facility class (functionally identical to x)
        model.z = pyo.Var(range(num_facilities), range(num_sites), domain=pyo.Binary, initialize=init_x_wrapper)
    model.lat = pyo.Var(range(num_facilities), domain=pyo.Reals, bounds=(lat_min, lat_max), initialize=random_lat_init)
    model.lon = pyo.Var(range(num_facilities), domain=pyo.Reals, bounds=(lon_min, lon_max), initialize=random_lon_init)

    model.assignment_constraint = pyo.ConstraintList()
    for j in range(num_sites):
        if is_hybrid:
            model.assignment_constraint.add(sum(model.x[i, j] for i in range(num_facilities)) + model.y[j] <= 1)
            model.assignment_constraint.add(sum(model.z[i, j] for i in range(num_facilities)) <= model.y[j])
        else:
            model.assignment_constraint.add(sum(model.x[i, j] for i in range(num_facilities)) <= 1)


    model.coverage_constraint = pyo.Constraint(rule=coverage_rule)

    return model

def model_init_site_locs_only(num_facilities, num_sites):

    model = pyo.ConcreteModel()


    def init_x_wrapper(model, i, j):
        return modeling_functions.initialize_x(model, i, j, num_facilities=num_facilities, num_sites=num_sites)

    def init_z_wrapper(model, i, j):
        return modeling_functions.initialize_z(model, i, j, num_facilities=num_facilities, num_sites=num_sites)

    model.x = pyo.Var(range(num_facilities), range(num_sites), domain=pyo.Binary, initialize=init_x_wrapper)
    model.z = pyo.Var(range(num_facilities), range(num_sites), domain=pyo.Binary, initialize=init_z_wrapper)

    model.assignment_constraint = pyo.ConstraintList()
    for j in range(num_sites):
        model.assignment_constraint.add(sum(model.x[i, j] for i in range(num_facilities)) <= 1)

    def facility_assignment_rule(model, i):
        return sum(model.z[i, j] for j in range(num_sites)) == 1

    model.facility_assignment = pyo.Constraint(range(num_facilities), rule=facility_assignment_rule)
    model.coverage_constraint = pyo.Constraint(
        expr=sum(model.x[i, j] for i in range(num_facilities) for j in range(num_sites)) >= 1 * num_sites)

    return model

