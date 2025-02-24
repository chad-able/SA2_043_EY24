import pyomo.environ as pyo
import logging
from pyomo.util.infeasible import log_infeasible_constraints
import random
import math

def inverse_order_of_magnitude(number):
    if number == 0:
        return "undefined"  # The order of magnitude for zero is undefined
    magnitude = math.floor(math.log10(abs(number)))
    return 10**(-magnitude)

def solver(model_type, executable_path):
    model_types = ['glpk', 'ipopt', 'bonmin', 'couenne', 'mindtpy']
    if model_type in model_types:
        if executable_path:
            solver = pyo.SolverFactory(model_type, executable=executable_path)
            return solver
        else:
            solver = pyo.SolverFactory(model_type)
    else:
        print("Model type not in valid model types")

def mindtpy_solve(mip_model_type, mip_executable_path, nlp_model_type, nlp_executable_path):
    return True

def log_pyomo_infeasible_constraints(model_instance):
    # Create a logger object with DEBUG level
    logging_logger = logging.getLogger()
    logging_logger.setLevel(logging.DEBUG)
    # Create a console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # add the handler to the logger
    logging_logger.addHandler(ch)
    # Log the infeasible constraints of pyomo object
    print("Displaying Infeasible Constraints")
    log_infeasible_constraints(model_instance, log_expression=True,
                         log_variables=True, logger=logging_logger)


def initialize_x(model, i, j, num_facilities, num_sites, is_hybrid):
    # Ensure at most one facility is assigned to each site
    # Pre-compute assignments for each site
    if is_hybrid:
        if not hasattr(model, 'x_assignments'):
            model.x_assignments = {}
            for j in range(num_sites):
                if model.y_assignments[j] == 0:
                    model.x_assignments[j] = random.randint(0, num_facilities - 1)
                else:
                    model.x_assignments[j] = None
        return 1 if model.x_assignments[j] == i else 0
    else:
        if not hasattr(model, 'x_assignments'):
            # Randomly assign one facility to each site
            model.x_assignments = {j: random.randint(0, num_facilities - 1) for j in range(num_sites)}

    # Assign 1 to the selected facility, 0 to others
    return 1 if model.x_assignments[j] == i else 0



def initialize_y(model, j, num_sites):
    # Precompute assignments only once and store them as an attribute of the model.
    if not hasattr(model, 'y_assignments'):
        # For each site j, randomly choose between 0 (centralized) and 1 (on-site)
        model.y_assignments = {j: random.choice([0, 1]) for j in range(num_sites)}
    return model.y_assignments[j]

def initialize_y_debug(model, j, num_sites):
    # Precompute assignments only once and store them as an attribute of the model.
    if not hasattr(model, 'y_assignments'):
        # For each site j, randomly choose between 0 (centralized) and 1 (on-site)
        model.y_assignments = {j: 1 for j in range(num_sites)}
    return model.y_assignments[j]
def initialize_z(model, i, j, num_facilities, num_sites, is_hybrid):

    if is_hybrid:
        if not hasattr(model, 'z_assignments'):
            # For each site that uses on-site treatment, randomly choose a solids facility.
            model.z_assignments = {}
            for j in range(num_sites):
                # Only assign if the precomputed y for this site is 1 (on-site)
                if model.y_assignments[j] == 1:
                    model.z_assignments[j] = random.randint(0, num_facilities - 1)
                else:
                    model.z_assignments[j] = None
        return 1 if model.z_assignments[j] == i else 0

    else:
        if not hasattr(model, 'z_assignments'):
            # Randomly assign each facility to a site
            model.z_assignments = {i: random.randint(0, num_sites - 1) for i in range(num_facilities)}

        # Assign 1 to the selected site, 0 to others
        return 1 if model.z_assignments[i] == j else 0

def logistic(x, n):
    log = x**n/(1+x**n)
    return log

def sigmoid_shifted(x, n, t):
    sig = 1/(1+pyo.exp(-n*(x-t)))
    return sig



