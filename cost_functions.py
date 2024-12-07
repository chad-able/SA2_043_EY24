import numpy as np
from distance_functions import haversine
import modeling_functions

def transportation_cost(dist,flow_rate):
    #assume constant trucking, 60 mph
    truck_capacity = 110*42 #initially in barrels, converted to gallons, from Pareto (p_delta_Truck in strategic_produced_water_optimization)
    truck_hourly_rate = 95 #$/hr, from Pareto, rates vary from 90 to 110
    truck_speed = 60*1609.34 # starts in mph, convert to meters per hour
    truck_time = 6*modeling_functions.sigmoid_shifted(dist, 10, 1)+dist/truck_speed # assume 6 hour minimum trucking time for loading/unloading, include travel time
    # truck_time = dist/truck_speed
    truck_num = np.ceil(flow_rate/truck_capacity*60*24*365) # number of trucks required per year
    truck_cost_per_year = truck_hourly_rate*truck_time*truck_num
    return truck_cost_per_year

def treatment_capex(total_flow_rate):
    # generic, assumed 1000 $/bbl/day from PARETO's treatment technology matrix, returns just $
    return total_flow_rate*1000*60*24/42

def treatment_opex(total_flow_rate):
    # generic, assumed $1/bbl feed from PARETO's treatment technology matrix, returns $/year
    return total_flow_rate*1*60*24*365/42

def annualized_cost(t_capex,t_opex,trucking,ut,crf,total_flow_rate):
    # returns $/bbl
    # ut = utilization, crf = capital recovery factor
    cost_annual = (crf*t_capex+t_opex+trucking)/(ut*total_flow_rate*60*24*365/42)
    return cost_annual

def cost_single_facility_any_annual(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0

    for j in range(num_sites):
        x_ij = model.x[index, j]
        t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex(total_flow_rate)
    opex = treatment_opex(total_flow_rate)
    annual_cost = annualized_cost(capex, opex, trans_cost, 1, 0.0769, total_flow_rate)
    # print('annual cost is ', annual_cost)

    return annual_cost

def cost_single_facility_any(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0

    for j in range(num_sites):
        x_ij = model.x[index, j]
        t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex(total_flow_rate)
    opex = treatment_opex(total_flow_rate)

    return capex, opex, trans_cost, total_flow_rate

def cost_single_facility_site_annual(model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0
    lat, lon = facility_coordinates(model, index, site_coordinates, num_sites)

    for j in range(num_sites):
        z_ij = model.z[index, j]
        x_ij = model.x[index, j]
        # t_cost = - (z_ij - 1) * x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        t_cost = x_ij * transportation_cost(
            haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex(total_flow_rate)
    opex = treatment_opex(total_flow_rate)
    annual_cost = annualized_cost(capex, opex, trans_cost, 1, 0.0769, total_flow_rate)
    # print('annual cost is ', annual_cost)

    return annual_cost

def cost_single_facility_site(model, index, num_sites, site_coordinates, flow_rate_data, h_approx):
    trans_cost = 0
    total_flow_rate = 0
    lat, lon = facility_coordinates(model, index, site_coordinates, num_sites)

    for j in range(num_sites):
        z_ij = model.z[index, j]
        x_ij = model.x[index, j]
        t_cost = - (z_ij - 1) * x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex(total_flow_rate)
    opex = treatment_opex(total_flow_rate)
    #annual_cost = annualized_cost(capex, opex, trans_cost, 1, 0.0769, total_flow_rate)
    # print('annual cost is ', annual_cost)

    return capex, opex, trans_cost, total_flow_rate

def facility_obj_per_annual(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag):
    total_cost = 0
    if not sites_flag:
        for i in range(num_facilities):
            # print('facility number ', i)
            cost_fac = cost_single_facility_any(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
            total_cost += cost_fac

    else:
        for i in range(num_facilities):
            cost_fac = cost_single_facility_site(model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
            total_cost += cost_fac

    # print('total_cost is', total_cost)
    return total_cost

def facility_obj(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag):
    total_cap = 0
    total_op = 0
    total_truck = 0
    total_flow = 0
    if not sites_flag:
        for i in range(num_facilities):
            # print('facility number ', i)
            capex, opex, trucking, flow_rate = cost_single_facility_any(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
            total_cap += capex
            total_op += opex
            total_truck += trucking
            total_flow += flow_rate

    else:
        for i in range(num_facilities):
            capex, opex, trucking, flow_rate = cost_single_facility_site(model, i, num_sites, site_coordinates, flow_rate_data, h_approx)
            total_cap += capex
            total_op += opex
            total_truck += trucking
            total_flow += flow_rate

    # print('total_cost is', total_cost)
    annual_cost = annualized_cost(total_cap, total_op, total_truck, 1, 0.0769, total_flow)
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