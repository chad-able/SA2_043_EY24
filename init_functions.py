import pandas as pd
import geopandas as gpd
import random
import pyomo.environ as pyo
import modeling_functions
import re


# TODO
# add functionality to separate out single play, list of plays
# drop wells that do not meet certain flow rate requirements
# Add in NEWTS concentration data, match to shale play
# Fix shale play relations

def read_well_data(filename, **gpd_args):
    layer = 1
    rows = 100
    if 'layer' in gpd_args:
        layer = gpd_args['layer']
    if 'rows' in gpd_args:
        rows = gpd_args['rows']
    gpdf = gpd.read_file(filename, layer=layer, rows=rows)
    return gpdf

def read_conc_data(filename):
    df = pd.read_csv(filename, encoding='ISO-8859-1', low_memory = False)
    return df

def filter_conc_data(df, plays):
    pattern = '|'.join(re.escape(play) for play in plays)

    df = df[df[['BASIN', 'FIELD', 'FORMATION']].apply(
        lambda row: row.str.contains(pattern, case=False, na=False).any(), axis=1
    )]

    return df


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

    water_locations_with_shales.to_csv("water_locations.csv", index=False)

    return restored_gpdf

def flow_rate_dump(gpdf):
    return gpdf

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


def model_init_any_location(lat_min, lat_max, lon_min, lon_max, num_facilities, num_sites):
    model = pyo.ConcreteModel()

    def random_lat_init(model, i):
        return random.uniform(lat_min, lat_max)

    def random_lon_init(model, i):
        return random.uniform(lon_min, lon_max)

    def init_x_wrapper(model, i, j):
        return modeling_functions.initialize_x(model, i, j, num_facilities=num_facilities, num_sites=num_sites)

    def coverage_rule(model):
        return sum(model.x[i, j] for i in range(num_facilities) for j in range(num_sites)) >= num_sites

    model.x = pyo.Var(range(num_facilities), range(num_sites), domain=pyo.Binary, initialize=init_x_wrapper)
    model.lat = pyo.Var(range(num_facilities), domain=pyo.Reals, bounds=(lat_min, lat_max), initialize=random_lat_init)
    model.lon = pyo.Var(range(num_facilities), domain=pyo.Reals, bounds=(lon_min, lon_max), initialize=random_lon_init)

    model.assignment_constraint = pyo.ConstraintList()
    for j in range(num_sites):
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