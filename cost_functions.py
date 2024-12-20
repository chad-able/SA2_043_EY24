import numpy as np
from distance_functions import haversine, linear_distance, euclidean_distance
import modeling_functions
# TODO
# separate cost functions for specific shale play (may require separate files)
# Modular + centralized, solids and liquid transport
# Modular relationships for volume reduction?
def transportation_cost(dist,flow_rate):
    # assume constant trucking, 60 mph
    # initially in barrels, converted to gallons, from Pareto (p_delta_Truck in strategic_produced_water_optimization, 110 from build_utils)
    truck_capacity = 110*42
    # $/hr, from Pareto, rates vary from 90 to 110 (strategic_toy_case_study, tab TruckingHourlyCost)
    truck_hourly_rate = 95
    # starts in mph, convert to meters per hour
    truck_speed = 60*1609.34
    # assume 6 hour minimum trucking time for loading/unloading, include travel time
    # sigmoid function used to bring truck_time to 0 at distance = 0
    truck_time = 6*modeling_functions.sigmoid_shifted(dist, 30, 1)+dist/truck_speed
    # number of trucks required per year
    truck_num = np.ceil(flow_rate/truck_capacity*60*24*365)
    truck_cost_per_year = truck_hourly_rate*truck_time*truck_num
    return truck_cost_per_year

def transportation_cost_solids(dist,flow_rate,density,concentration):
    # density in kg/m3
    # concentration in mg/L

    # mass of solute, in milligrams per minute
    mass_solute = concentration * flow_rate * 3.78541

    # volume of solute in gallons per minute
    vol_solute = mass_solute / 1000000 / density * 264.172

    # compute based on transporting equivalent volume of liquid
    trans_cost = transportation_cost(dist, vol_solute)

    return trans_cost




def treatment_capex_central(total_flow_rate):
    # generic, assumed 1000 $/bbl/day from PARETO's treatment technology matrix, returns just $
    return total_flow_rate*1000*60*24/42

def treatment_capex_central_scaled(total_flow_rate):
    # base capacity of the central plant is most commonly 5079 bbl/day from PARETO's treatment technology matrix
    # attempts to scale costs vs. this base capacity

    initial_cost = treatment_capex_central(total_flow_rate)
    scaled_cost = initial_cost*((total_flow_rate*60*24/42)/5079)**0.6
    return scaled_cost

def treatment_capex_modular(total_flow_rate):
    # generic, assumed 3 times higher than centralized costs, returns just $
    return total_flow_rate*3*1000*60*24/42


def treatment_opex_central(total_flow_rate):
    # generic, assumed $1/bbl feed from PARETO's treatment technology matrix, returns $/year
    return total_flow_rate*1*60*24*365/42

def treatment_opex_central_scaled(total_flow_rate):
    # base capacity of the central plant is most commonly 5079 bbl/day from PARETO's treatment technology matrix
    # attempts to scale costs vs. this base capacity

    initial_cost = treatment_opex_central(total_flow_rate)
    scaled_cost = initial_cost*((total_flow_rate*60*24/42)/5079)**0.6
    return scaled_cost

def treatment_opex_modular(total_flow_rate):
    # generic, assumed 5 times higher than centralized costs, returns $/year
    return total_flow_rate*5*60*24*365/42

def annualized_cost(t_capex,t_opex,trucking,total_flow_rate,ut=1,crf=0.0769):
    # returns $/bbl
    # ut = utilization, crf = capital recovery factor
    cost_annual = (crf*t_capex+t_opex+trucking)/(ut*total_flow_rate*60*24*365/42)
    return cost_annual


def cost_single_facility_any_central(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0

    for j in range(num_sites):
        x_ij = model.x[index, j]
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        t_cost = x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex_central_scaled(total_flow_rate)
    opex = treatment_opex_central_scaled(total_flow_rate)

    return capex, opex, trans_cost, total_flow_rate

def cost_sites_modular(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx, density, concentration):
    trans_cost = 0
    total_flow_rate = 0
    for j in range(num_sites):
        x_ij = model.x[index, j]
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        t_cost = x_ij * transportation_cost_solids(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0], density, concentration)
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow
    capex = treatment_capex_modular(total_flow_rate)
    opex = treatment_opex_modular(total_flow_rate)

    return capex, opex, trans_cost, total_flow_rate

def annual_cost_modular(num_sites, flow_rate_data):
    total_flow = 0
    total_cap = 0
    total_op = 0
    total_truck = 0
    for j in range(num_sites):
        capex, opex, trans_cost = cost_sites_modular(flow_rate_data[j][0])
        total_cap += capex
        total_op += opex
        total_truck += trans_cost
        total_flow += flow_rate_data[j][0]

    cost_annual = annualized_cost(total_cap, total_op, total_truck, total_flow, 1, 0.0769)
    return cost_annual

def cost_single_facility_site_central(model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0
    lat, lon = facility_coordinates(model, index, site_coordinates, num_sites)

    for j in range(num_sites):
        z_ij = model.z[index, j]
        x_ij = model.x[index, j]
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        t_cost = x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex_central_scaled(total_flow_rate)
    opex = treatment_opex_central_scaled(total_flow_rate)

    return capex, opex, trans_cost, total_flow_rate


def facility_obj(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag, mod_flag, density = 534, concentration = 10):
    total_cap = 0
    total_op = 0
    total_truck = 0
    total_flow = 0
    if mod_flag:
        for i in range(num_facilities):
            capex, opex, trucking, flow_rate = cost_sites_modular(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx, density, concentration)
            total_cap += capex
            total_op += opex
            total_truck += trucking
            total_flow += flow_rate
    else:
        if not sites_flag:
            for i in range(num_facilities):
                # print('facility number ', i)
                capex, opex, trucking, flow_rate = cost_single_facility_any_central(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
                total_cap += capex
                total_op += opex
                total_truck += trucking
                total_flow += flow_rate

        else:
            for i in range(num_facilities):
                capex, opex, trucking, flow_rate = cost_single_facility_site_central(model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
                total_cap += capex
                total_op += opex
                total_truck += trucking
                total_flow += flow_rate

    # print('total_cost is', total_cost)
    annual_cost = annualized_cost(total_cap, total_op, total_truck, total_flow, 1, 0.0769)
    return annual_cost

def transportation_cost_total(model):
    return sum((model.x[i, j] * transportation_cost(haversine(model.lat[i], model.lon[i], site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0]) for i in range(num_facilities) for j in range(num_sites)))


def x_fac(model, i, num_sites):
    return [model.x[i,j] for j in range(num_sites)]

def facility_coordinates(model, i, site_coordinates, num_sites):
    x_coord = sum(site_coordinates[j][0] * model.z[i, j] for j in range(num_sites))
    y_coord = sum(site_coordinates[j][1] * model.z[i, j] for j in range(num_sites))
    return x_coord, y_coord

def facility_coordinates_vals(model, i, site_coordinates, num_sites):
    x_coord = sum(site_coordinates[j][0] * model.z[i, j].value for j in range(num_sites))
    y_coord = sum(site_coordinates[j][1] * model.z[i, j].value for j in range(num_sites))
    return x_coord, y_coord