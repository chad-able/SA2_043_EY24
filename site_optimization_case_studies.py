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
import re



def main():
    facility_array = [2]
    attempts_per = 1
    random_seed_use = False
    random_seed = 43424267667

    if random_seed_use:
        random.seed(random_seed)




    filename = "E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb"

    gpd_args = {
        'layer': 1,
        'rows': 1000
    }

    # initial data read
    huc_data = init_functions.read_data(filename, **gpd_args)

    #addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
    huc_data = init_functions.centroids(huc_data, 0)






    huc_data = init_functions.flow_rate(huc_data)

    huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)

    print(huc_data.head())

    # shale_filename = "E:/codes/RC24/PWMapping/SedimentaryBasins_US_EIA/Lower_48_Sedimentary_Basins.shp"
    shale_filename = "E:/codes/RC24/PWMapping/TightOil_ShaleGas_Plays_Lower48_EIA/TightOil_ShaleGas_Plays_Lower48_EIA.shp"
    huc_data = init_functions.shale_plays(huc_data, shale_filename)
    plays = ['Bakken', 'Eagle Ford', 'Permian']
    pattern = '|'.join(re.escape(play) for play in plays)
    # print(huc_data["Shale_play"])

    huc_data = huc_data[huc_data['Shale_play'].str.contains(pattern, case=False, na=False)]
    num_sites = len(huc_data)
    print('num_sites is ', num_sites)
    lat_max = huc_data['lat'].max()
    lat_min = huc_data['lat'].min()
    lon_max = huc_data['lon'].max()
    lon_min = huc_data['lon'].min()

    h_approx = True
    site_coordinates = list(zip(huc_data['lat'], huc_data['lon']))
    flow_rate_data = list(zip(huc_data["2022_flow_gpm"]))
    # if sites_flag = True, use only sites as locations
    # Otherwise, use any valid location

    for fnum in range(len(facility_array)):
        num_facilities = facility_array[fnum]
        sites_flag = False
        mod_flag = False

        for run in range(attempts_per):
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


            model.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
            model.scaling_factor[model.objective] = modeling_functions.inverse_order_of_magnitude(testval)

            if not sites_flag:
                model.scaling_factor[model.lat] = modeling_functions.inverse_order_of_magnitude(lat_max)
                model.scaling_factor[model.lon] = modeling_functions.inverse_order_of_magnitude(lon_max)

            scaled_model = pyo.TransformationFactory('core.scale_model').create_using(model)

            solvers = ['couenne', 'bonmin', 'glpk']

            selected_solver = solvers[0]

            executable_paths = {
                'couenne': "C:/Users/chada/AppData/Local/idaes/bin/couenne",
                'bonmin': "C:/Users/chada/AppData/Local/idaes/bin/bonmin",
                'glpk': 'E:\codes\glpk-4.65\w64\glpsol',
            }

            solver = modeling_functions.solver(selected_solver, executable_paths[selected_solver])


            if selected_solver == 'couenne':
                couenne_dict = {
                    'max_cpu_time': 600,
                    'max_iter': 1000,
                    'tol': 1e-3,
                    'bonmin.time_limit': 100,
                    'bonmin.node_limit': 10000,
                    'bonmin.resolve_on_small_infeasibility': 1,
                    'bonmin.cutoff': 10,
                    'art_lower': 0,
                    'art_cutoff': 10,
                    'feas_tolerance': 1e-4,
                    'bonmin.algorithm': "B-OA",
                    'bonmin.node_comparison': "dynamic",
                    'bonmin.integer_tolerance': 1e-3,
                    'bonmin.num_retry_unsolved_random_point': 5,
                    'bonmin.acceptable_tol': 1e-3
                }

                with open('couenne.opt', 'w') as file:
                    for key, value in couenne_dict.items():
                        write_str = f"{key} {value}"
                        file.write(write_str + "\n")

                result = solver.solve(scaled_model, tee=True)


            if selected_solver == 'bonmin':


                result = solver.solve(scaled_model, tee=True, options={
                    'max_iter': 1000,
                    'bonmin.time_limit': 600,
                    'tol': 1e-4,
                    'acceptable_tol': 1e-3,
                    'bonmin.algorithm': 'B-BB',
                    'bonmin.oa_decomposition': 'yes',
                    'bonmin.node_comparison': 'dynamic',
                    'bonmin.num_retry_unsolved_random_point': 5,
                    'bonmin.random_generator_seed': -1
                })


            modeling_functions.log_pyomo_infeasible_constraints(scaled_model)

            pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_model, model)


            print_file = 'facility_output' + str(num_facilities) + '_' + str(run) + '.csv'
            newval = pyo.value(objective_f(model))
            finalprint = print_functions.print_final(print_array, newval, model, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag, print_file)


if __name__ == "__main__":
    m = main()



