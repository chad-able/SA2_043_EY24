import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry import Polygon
from scipy import optimize
import numpy as np
import pyomo.environ as pyo
import math
import random
from pyomo.util.infeasible import log_infeasible_constraints
import logging
import init_functions
import distance_functions
import cost_functions
import modeling_functions
from functools import partial
import csv
import print_functions




num_facilities = 2
random_seed = 443

#random.seed(random_seed)
# huc8 = pd.read_csv("huc8_summary.csv")
# for i in range(len(huc8.columns)):
#     print(huc8.columns[i])

filename = ("E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb")

gpd_args = {
    'layer': 1,
    'rows': 250
}

# initial data read
huc_data = init_functions.read_data(filename, **gpd_args)

#addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
huc_data = init_functions.centroids(huc_data, 0)


lat_max = huc_data['lat'].max()
lat_min = huc_data['lat'].min()
lon_max = huc_data['lon'].max()
lon_min = huc_data['lon'].min()

h_approx = True

site_coordinates = list(zip(huc_data['lat'], huc_data['lon']))

huc_data = init_functions.flow_rate(huc_data)

huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)
flow_rate_data = list(zip(huc_data["2022_flow_gpm"]))
print(huc_data.head())
num_sites = len(huc_data)
shale_filename = "E:/codes/RC24/PWMapping/SedimentaryBasins_US_EIA/Lower_48_Sedimentary_Basins.shp"
#huc_data = init_functions.shale_plays(huc_data, shale_filename)

#print(huc_data["Shale_play"])

# if sites_flag = True, use only sites as locations
# Otherwise, use any valid location


sites_flag = False
mod_flag = False

if not sites_flag:
    model = init_functions.model_init_any_location(lat_min, lat_max, lon_min, lon_max, num_facilities, num_sites)

else:
    model = init_functions.model_init_site_locs_only(num_facilities, num_sites)




def objective_f(model):
    return cost_functions.facility_obj(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag, mod_flag)

print(flow_rate_data[0][0])
# print("Modular cost is ", cost_functions.annual_cost_modular(num_sites, flow_rate_data))
testval = pyo.value(objective_f(model))
print_array = print_functions.print_initial(testval, model, num_facilities, num_sites, site_coordinates, sites_flag)


model.objective = pyo.Objective(rule=objective_f, sense=pyo.minimize)
#model.pprint()

model.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
model.scaling_factor[model.objective] = modeling_functions.inverse_order_of_magnitude(testval)

if not sites_flag:
    model.scaling_factor[model.lat] = modeling_functions.inverse_order_of_magnitude(lat_max)
    model.scaling_factor[model.lon] = modeling_functions.inverse_order_of_magnitude(lon_max)

scaled_model = pyo.TransformationFactory('core.scale_model').create_using(model)


executable_paths = {
    'couenne': "C:/Users/chada/AppData/Local/idaes/bin/couenne",
    'bonmin': "C:/Users/chada/AppData/Local/idaes/bin/bonmin",
    'glpk': 'E:\codes\glpk-4.65\w64\glpsol',
}
#nlp_solver = pyo.SolverFactory('ipopt')
#solver = pyo.SolverFactory('cbc')
solver = modeling_functions.solver('couenne', executable_paths['couenne'])

couenne_dict = {
    'max_cpu_time': 600,
    'max_iter': 1000,
    'tol': 1e-4,
    'bonmin.time_limit': 100,
    'bonmin.node_limit': 1000,
    'bonmin.resolve_on_small_infeasibility': 1,
    'bonmin.cutoff': 10,
    'art_lower': 0,
    'art_cutoff': 10,
    'bonmin.algorithm': "B-Hyb",
    'bonmin.node_comparison': "dynamic",
    'bonmin.integer_tolerance': 1e-3,
    'bonmin.acceptable_tol': 1e-3
}

with open('couenne.opt', 'w') as file:
    for key, value in couenne_dict.items():
        write_str = f"{key} {value}"
        file.write(write_str + "\n")




result = solver.solve(scaled_model, tee=True, options={
    'max_iter': 1000,
    'bonmin.time_limit': 600,
    'tol': 1e-4,
    'acceptable_tol': 1e-3,
})

# result = solver.solve(scaled_model, tee=True, options={
#     'max_iter': 100,
#     'max_cpu_time': 200,
#     'tol': 1e-5,
#     'acceptable_tol': 1e-3,
#     'bonmin.time_limit': 10,
#     'bonmin.node_limit': 300,
# })

# result = solver.solve(scaled_model, tee=True)


modeling_functions.log_pyomo_infeasible_constraints(scaled_model)

pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_model, model)

#result = solver.solve(model, mip_solver = 'cbc', nlp_solver = 'ipopt', tee=True, mip_solver_tee=True, nlp_solver_tee=True)
#

newval = pyo.value(objective_f(model))
finalprint = print_functions.print_final(print_array, newval, model, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag)






# solver = pyo.SolverFactory('ipopt')
# result = solver.solve(model, tee=True)
# def sample_treatment_distance(p):
# #create sample treatment facility
#     treatment_test = Point(p[0],p[1])
#     huc_data["2022_flow_gpm"] = huc_data["F2022_WaterProd_BBL_sum_1"].apply(flow_rate_calc)
#     total_flow_rate = huc_data["2022_flow_gpm"].sum()
#     huc_data["treatment_distance"] = huc_data.apply(lambda row: row.geometry_centroids.distance(treatment_test), axis=1)
#     huc_data["trucking_cost"] = huc_data.apply(lambda row: transportation_cost(row["treatment_distance"], row["2022_flow_gpm"]), axis=1)
#     trucking = huc_data["trucking_cost"].sum()
#     capex = treatment_capex(total_flow_rate)
#     opex = treatment_opex(total_flow_rate)
#     annual_cost = annualized_cost(capex,opex,trucking,1,0.0769,total_flow_rate)

    #
    # return annual_cost

# init_points = (-1000000,1000000)
# min_point = optimize.minimize(sample_treatment_distance, init_points, method="Nelder-Mead")
# print(min_point)

# centroid_test1 = huc8_data["geometry_centroids"].to_crs('+proj=cea')
# centroid_test1 = huc8_data["geometry_centroids"]
# print(centroid_test1)
# poly = Polygon([[p.x, p.y] for p in centroid_test1])
# centroid_test2 = poly.centroid
# print(centroid_test2)

# plot = test.plot(column="F2022_WaterProd_BBL_sum_1")
# plt.show()