import numpy as np
from distance_functions import haversine, linear_distance, euclidean_distance
import modeling_functions
# TODO
# separate cost functions for specific shale play (may require separate files)
# Modular + centralized, solids and liquid transport
# Modular relationships for volume reduction?
def transportation_cost(dist,flow_rate, **params):
    if params['use_pareto']:
        # assume constant trucking
        # initially in barrels, converted to gallons
        truck_capacity = params['truck_capacity_liq']*42
        # $/hr
        truck_hourly_rate = params['truck_hourly_rate']
        # starts in mph, convert to km per hour
        truck_speed = params['truck_driving_speed']*1.60934
        # sigmoid function used to bring truck_time to 0 at distance = 0
        truck_time = params['truck_loading_time']*modeling_functions.sigmoid_shifted(dist, 30, 1)+dist/truck_speed
        # number of trucks required per year
        truck_num = np.ceil(flow_rate/truck_capacity*60*24*365)
        truck_cost_per_year = truck_hourly_rate*truck_time*truck_num
    else:
        # much simpler relation, carries in cost in $/bbl/mile, converts using flow rate and distance, generates $/year
        truck_cost_per_year = params['truck_variable_cost']*dist/1.60934*flow_rate/42*60*24*365
    return truck_cost_per_year

def transportation_cost_solids(dist,flow_rate,**params):
    # density in kg/m3
    # concentration in mg/L

    # mass of solute, in milligrams per minute
    mass_solute = params['concentration'] * flow_rate * 3.78541

    # volume of solute in gallons per minute
    vol_solute = mass_solute / 1000000 / params['density'] * 264.172

    # compute based on transporting equivalent volume of liquid
    trans_cost = transportation_cost(dist, vol_solute, **params['trucking'])

    return trans_cost

def solids_conversions(flow_rate, **params):
    # density in kg/m3
    # concentration in mg/L

    # mass of solute, in milligrams per minute

    mass_solute = params['concentration'] * flow_rate * 3.78541

    # volume of solute in gallons per minute

    vol_solute = mass_solute / 1000000 / params['density'] * 264.172

    return mass_solute, vol_solute



def treatment_capex_central(total_flow_rate, **params):
    # central cap cost in $/bbl/day, returns $
    return total_flow_rate*params['central_cap_cost_init']*60*24/42

def treatment_capex_central_scaled(total_flow_rate, **params):
    # attempts to scale costs vs. base capacity, simple exponential curve, returns $
    eps = 1e-6  # to ensure that power function is never 0

    initial_cost = treatment_capex_central(total_flow_rate, **params)
    scaled_cost = initial_cost*(((total_flow_rate+eps)*60*24/42)/params['base_central_capacity'])**params['central_cap_scale_exp']
    return scaled_cost

def treatment_capex_modular(total_flow_rate, **params):
    # modular cap cost in $/bbl/day, returns $
    return total_flow_rate*params['modular_cap_cost_init']*60*24/42


def treatment_opex_central(total_flow_rate, **params):
    # central op cost in $/bbl, returns $/year
    return total_flow_rate*params['central_op_cost_init']*60*24*365/42

def treatment_opex_central_scaled(total_flow_rate, **params):
    # attempts to scale costs vs. base capacity, simple exponential curve, returns $/year
    eps = 1e-6  # to ensure that power function is never 0

    initial_cost = treatment_opex_central(total_flow_rate, **params)
    scaled_cost = initial_cost*(((total_flow_rate+eps)*60*24/42)/params['base_central_capacity'])**params['central_op_scale_exp']
    return scaled_cost

def treatment_capex_solids(total_solids, **params):
    # $/tonne/day? total solids in mg/min. Returns $
    return total_solids * params['solids_cap_cost_init'] * 60 * 24 / 1000000000

def treatment_capex_solids_scaled(total_solids, **params):
    # attempts to scale costs vs. base capacity, simple exponential curve, returns $
    eps = 1e-6  # to ensure that power function is never 0
    initial_cost = treatment_capex_solids(total_solids, **params)
    scaled_cost = initial_cost*(((total_solids+eps)*60*24/1000000000)/params['base_solids_capacity'])**params['solids_cap_scale_exp']
    return scaled_cost

def treatment_opex_solids(total_solids, **params):
    # central op cost in $/bbl, returns $/year
    return total_solids*params['solids_op_cost_init']*60*24*365/1000000000

def treatment_opex_solids_scaled(total_solids, **params):
    # attempts to scale costs vs. base capacity, simple exponential curve, returns $/year
    eps = 1e-6  # to ensure that power function is never 0

    initial_cost = treatment_opex_solids(total_solids, **params)
    scaled_cost = initial_cost*(((total_solids+eps)*60*24/1000000000)/params['base_solids_capacity'])**params['solids_op_scale_exp']
    return scaled_cost

def treatment_opex_modular(total_flow_rate, **params):
    # generic, assumed 5 times higher than centralized costs, returns $/year
    return total_flow_rate*params['modular_op_cost_init']*60*24*365/42

def annualized_cost(t_capex,t_opex,trucking,total_flow_rate,ut=1,crf=0.0769):
    # returns $/bbl
    cost_annual = (crf*t_capex+t_opex+trucking)/(ut*total_flow_rate*60*24*365/42)
    return cost_annual


def cost_single_facility_any_central(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params):
    trans_cost = 0
    total_flow_rate = 0
    truck_params = params['trucking']
    for j in range(num_sites):
        x_ij = model.x[index, j]
        if is_hybrid:
            y_ij = model.y[j]
        else:
            y_ij = 0
        # t_cost = (1 - y_ij) * x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])

        t_cost = (1 - y_ij) * x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0], **truck_params)
        trans_cost += t_cost
        flow = (1 - y_ij) * x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    cost_params = params['costing']
    capex = treatment_capex_central_scaled(total_flow_rate, **cost_params)
    opex = treatment_opex_central_scaled(total_flow_rate, **cost_params)

    return capex, opex, trans_cost, total_flow_rate

def cost_sites_solids(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params):
    trans_cost = 0
    total_flow_rate = 0
    total_mass = 0
    truck_params = params['trucking']
    for j in range(num_sites):

        if is_hybrid:
            y_ij = model.y[j]
            x_ij = model.z[index, j]
        else:
            y_ij = 1
            x_ij = model.x[index, j]
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        solids_mass, solids_flow = solids_conversions(flow_rate_data[j][0], **params)
        t_cost = y_ij * x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), solids_flow, **truck_params)
        trans_cost += t_cost
        flow = y_ij * x_ij * flow_rate_data[j][0]
        total_flow_rate += flow
        total_mass += y_ij * solids_mass
    costing_params = params['costing']
    capex = treatment_capex_solids_scaled(total_mass, **costing_params)
    opex = treatment_opex_solids_scaled(total_mass, **costing_params)

    return capex, opex, trans_cost, total_flow_rate

def total_cost_modular(num_sites, flow_rate_data, model, is_hybrid, **params):
    total_cap = 0
    total_op = 0
    cost_params = params['costing']
    if is_hybrid:
        for j in range(num_sites):
            y_ij = model.y[j]
            flow = (1 - y_ij) * flow_rate_data[j][0]
            total_cap += treatment_capex_modular(flow, **cost_params)
            total_op += treatment_opex_modular(flow, **cost_params)
    else:
        for j in range(num_sites):
            total_cap += treatment_capex_modular(flow_rate_data[j][0], **cost_params)
            total_op += treatment_opex_modular(flow_rate_data[j][0], **cost_params)

    return total_cap, total_op

def cost_single_facility_site_central(model, index, num_sites, site_coordinates, flow_rate_data, h_approx, **params):
    trans_cost = 0
    total_flow_rate = 0
    lat, lon = facility_coordinates(model, index, site_coordinates, num_sites)

    for j in range(num_sites):
        z_ij = model.z[index, j]
        x_ij = model.x[index, j]
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])
        t_cost = x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0], **params['trucking'])
        trans_cost += t_cost
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    capex = treatment_capex_central_scaled(total_flow_rate, **params['costing'])
    opex = treatment_opex_central_scaled(total_flow_rate, **params['costing'])

    return capex, opex, trans_cost, total_flow_rate

def cost_sites_hybrid(lat, lon, model, index, num_sites, site_coordinates, flow_rate_data, h_approx, **params):
    trans_cost = 0
    total_flow_rate = 0
    truck_params = params['trucking']
    for j in range(num_sites):
        x_ij = model.x[index, j]
        y_ij = model.y[index, j]
        z_ij = model.z[index, j]
        solids_mass, solids_flow = solids_conversions(flow_rate_data[j][0], **params)
        # t_cost = x_ij * transportation_cost(haversine(lat, lon, site_coordinates[j][0], site_coordinates[j][1], h_approx), flow_rate_data[j][0])

        t_cost_liq = (1 - y_ij) * x_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), flow_rate_data[j][0], **truck_params)
        t_cos_sol = y_ij * z_ij * transportation_cost(
            euclidean_distance(lat, lon, site_coordinates[j][0], site_coordinates[j][1]), solids_flow, **truck_params)
        trans_cost += t_cost_liq + t_cost_sol
        flow = x_ij * flow_rate_data[j][0]
        total_flow_rate += flow

    cost_params = params['costing']
    capex = treatment_capex_central_scaled(total_flow_rate, **cost_params)
    opex = treatment_opex_central_scaled(total_flow_rate, **cost_params)

    return capex, opex, trans_cost, total_flow_rate
def facility_obj(model, num_sites, num_facilities, site_coordinates, flow_rate_data, h_approx, sites_flag, mod_flag, is_hybrid, **params):
    total_cap = 0
    total_op = 0
    total_truck = 0
    total_flow = 0
    if is_hybrid:
        cap_modular, op_modular = total_cost_modular(num_sites, flow_rate_data, model, is_hybrid, **params)
        total_cap += cap_modular
        total_op += op_modular
        for i in range(num_facilities):
            liq_cap, liq_op, liq_truck, liq_flow = cost_single_facility_any_central(model.lat_c[i],model.lon_c[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params)
            sol_cap, sol_op, sol_truck, sol_flow = cost_sites_solids(model.lat_s[i],model.lon_s[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params)
            total_cap += liq_cap + sol_cap
            total_op += liq_op + sol_op
            total_truck += liq_truck + sol_truck
            total_flow += liq_flow + sol_flow
    else:
        if mod_flag:
            cap_modular, op_modular = total_cost_modular(num_sites, flow_rate_data, model, is_hybrid, **params)
            total_cap += cap_modular
            total_op += op_modular
            for i in range(num_facilities):
                capex, opex, trucking, flow_rate = cost_sites_solids(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params)
                total_cap += capex
                total_op += opex
                total_truck += trucking
                total_flow += flow_rate
        else:
            if not sites_flag:
                for i in range(num_facilities):
                    # print('facility number ', i)
                    capex, opex, trucking, flow_rate = cost_single_facility_any_central(model.lat[i],model.lon[i], model, i, num_sites, site_coordinates, flow_rate_data, h_approx, is_hybrid, **params)
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