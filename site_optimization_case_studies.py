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
from time import time



def main():
    allstart = time()
    facility_array = [2]
    attempts_per = 10
    random_seed_use = False
    random_seed = 43424267667

    if random_seed_use:
        random.seed(random_seed)




    well_filename = "E:/codes/RC24/PWMapping/NEWTS_Well_Summary_by_Hydrologic_Regions_and_Subbasins.gdb"

    gpd_args = {
        'layer': 1,
    }


    # initial data read
    huc_data = init_functions.read_well_data(well_filename, **gpd_args)

    #addition of centroids to dataframe, proj_flag: 0 corresponds to EPSG:4326
    huc_data = init_functions.centroids(huc_data, 0)






    huc_data = init_functions.flow_rate(huc_data)

    huc_data["2022_flow_gpm"].fillna(huc_data["2022_flow_gpm"].mean(), inplace=True)

    print(huc_data.head())

    shale_filename = "E:/codes/RC24/PWMapping/SedimentaryBasins_US_EIA/Lower_48_Sedimentary_Basins.shp"
    # shale_filename = "E:/codes/RC24/PWMapping/TightOil_ShaleGas_Plays_Lower48_EIA/TightOil_ShaleGas_Plays_Lower48_EIA.shp"


    newts_filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
    newts_data = init_functions.read_conc_data(newts_filename)
    huc_data = init_functions.combine_conc_data(huc_data, newts_data)

    plays = ['Permian']
    pattern = '|'.join(re.escape(play) for play in plays)
    # fix shale play relations
    huc_data = init_functions.shale_plays(huc_data, shale_filename)
    huc_data = init_functions.filter_by_shale_plays(huc_data, plays, replace_values=False)
    num_sites = len(huc_data)
    print('num_sites is ', num_sites)
    lat_max = huc_data['lat'].max()
    lat_min = huc_data['lat'].min()
    lon_max = huc_data['lon'].max()
    lon_min = huc_data['lon'].min()

    h_approx = True
    site_coordinates = list(zip(huc_data['lat'], huc_data['lon']))
    flow_rate_data = list(zip(huc_data["2022_flow_gpm"]))
    # need to eventually include array of li conc data
    # Li_conc = list(zip(huc_data['Li']))
    Li_conc = huc_data['Li'].median()
    # if sites_flag = True, use only sites as locations
    # Otherwise, use any valid location

    alt_trucking_params = {
        'use_pareto': False,
        'truck_variable_cost': 0.10 # $/bbl/mile, from https://doi.org/10.2118/157532-MS (range from 0.02 to 0.40, from 2012)
        # other sources for trucking cost include doi:10.1002/aic.14705, https://www.osti.gov/servlets/purl/1533719
    }

    trucking_params = {
        'truck_capacity_liq': 110, # in barrels, from Pareto (p_delta_Truck in strategic_produced_water_optimization, 110 from build_utils)
        'truck_hourly_rate': 95, # $/hr, from Pareto, rates vary from 90 to 110 (strategic_toy_case_study, tab TruckingHourlyCost)
        'truck_driving_speed': 60, # mph
        'truck_loading_time': 6 # assume 6 hour minimum trucking time for loading/unloading, include travel time
    }



    costing_params = {
        'central_cap_cost_init': 1000,  # $/bbl/day, assumed 1000 $/bbl/day from PARETO's treatment technology matrix
        'central_op_cost_init': 1,  # $/bbl, assumed $1/bbl feed from PARETO's treatment technology matrix
        'modular_cap_cost_init': 5000,  # $/bbl/day
        'modular_op_cost_init': 3,  # $/bbl
        'solids_cap_cost_init': 3000, # $/tonne/day?
        'solids_op_cost_init': 10, # $/tonne?
        'central_cap_scale_exp': 0.6,  # for scaling central capital costs
        'modular_cap_scale_exp': 0.8,  # for scaling modular capital costs (not currently supported)
        'solids_cap_scale_exp': 0.6, # for scaling solids treatment capital costs
        'central_op_scale_exp': 1,  # for scaling central operating costs
        'modular_op_scale_exp': 1,  # for scaling modular operating costs (not currently supported)
        'solids_op_scale_exp': 1, # for scaling solids treatment operating costs
        'base_central_capacity': 5079, # base capacity of the central plant is most commonly 5079 bbl/day from PARETO's treatment technology matrix
        'base_solids_capacity': 10,  # tonnes/day?
        'utilization': 1, # fraction of uptime, between 0 and 1
        'capital_recovery_factor': 0.0769 # capital recovery factor for plant life

    }

    texas_param_dict = {
        'density': 534, # density of solid Li
        'concentration': Li_conc, # mg/L Li in produced water, currently median of produced water values
        'costing': costing_params,
        'trucking': trucking_params # nested dict for trucking parameters

    }

    param_dict = texas_param_dict
    for fnum in range(len(facility_array)):
        num_facilities = facility_array[fnum]


        sites_flag = False
        mod_flag = False
        is_hybrid = True
        best_val = 1000000
        for run in range(attempts_per):
            if not sites_flag:
                model = init_functions.model_init_any_location(lat_min, lat_max, lon_min, lon_max, num_facilities, num_sites, is_hybrid)

            else:
                model = init_functions.model_init_site_locs_only(num_facilities, num_sites)




            def objective_f(model):
                return cost_functions.facility_obj(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag, mod_flag, is_hybrid, **param_dict)

            print(flow_rate_data[0][0])
            # print("Modular cost is ", cost_functions.annual_cost_modular(num_sites, flow_rate_data))
            testval = pyo.value(objective_f(model))
            print_array = print_functions.print_initial(testval, model, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag, is_hybrid)


            model.objective = pyo.Objective(rule=objective_f, sense=pyo.minimize)


            model.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
            model.scaling_factor[model.objective] = modeling_functions.inverse_order_of_magnitude(testval)

            if is_hybrid:
                model.scaling_factor[model.lat_c] = modeling_functions.inverse_order_of_magnitude(lat_max)
                model.scaling_factor[model.lon_c] = modeling_functions.inverse_order_of_magnitude(lon_max)
                model.scaling_factor[model.lat_s] = modeling_functions.inverse_order_of_magnitude(lat_max)
                model.scaling_factor[model.lon_s] = modeling_functions.inverse_order_of_magnitude(lon_max)
            elif not sites_flag:
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
            model.display()

            if selected_solver == 'couenne':
                couenne_dict = {
                    'acceptable_tol': 1e-3,
                    'max_cpu_time': 1800,
                    'max_iter': 10000,
                    'tol': 1e-2,
                    'bonmin.time_limit': 1000,
                    'bonmin.node_limit': 10000,
                    'bonmin.resolve_on_small_infeasibility': 3,
                    'bonmin.cutoff': testval,
                    'art_lower': 0,
                    'art_cutoff': testval,
                    'log_num_obbt_per_level': 20,
                    'log_num_local_optimization_per_level': 10,
                    'feas_tolerance': 1e-3,
                    'bonmin.algorithm': "B-OA",
                    'bonmin.node_comparison': "dynamic",
                    'bonmin.integer_tolerance': 1e-3,
                    'bonmin.num_retry_unsolved_random_point': 10,
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
            timestr = 'Total time: {}s'.format(round(time() - allstart), 3)
            print_debug = print_functions.print_final(timestr, print_array, newval, model, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag, is_hybrid, print_file)
            if newval < best_val:
                best_val = newval
                model_save = model
                print_save = print_array
                best_file = 'facility_output' + str(num_facilities) + '_' + str(run) + 'best' + '.csv'

        timestr = 'Total time: {}s'.format(round(time()-allstart),3)
        finalprint = print_functions.print_final(timestr, print_save, best_val, model_save, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag, is_hybrid, best_file)


if __name__ == "__main__":
    m = main()



