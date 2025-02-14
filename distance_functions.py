import pyomo.environ as pyo
import math




def asin_taylor(x, n_terms=5):
    """Only valid if x is in the range [-1, 1]."""

    # if x > 1:
    #     raise ValueError("x must be in the range [-1, 1].")
    #
    # if x < -1:
    #     raise ValueError("x must be in the range [-1, 1].")

    result = 0
    for n in range(n_terms):
        coefficient = math.factorial(2 * n) / (4**n * (math.factorial(n)**2) * (2 * n + 1))
        result += coefficient * (x**(2 * n + 1))
    return result

def haversine(lat1, lon1, lat2, lon2, approx):
    # Convert degrees to radians
    lat1 = lat1 * math.pi / 180
    lon1 = lon1 * math.pi / 180
    lat2 = lat2 * math.pi / 180
    lon2 = lon2 * math.pi / 180

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    eps = 1e-6 #to ensure that sqrt is never 0
    a = pyo.sin(dlat / 2) ** 2 + pyo.cos(lat1) * pyo.cos(lat2) * pyo.sin(dlon / 2) ** 2
    if approx:
        c = 2 * asin_taylor(pyo.sqrt(a + eps))
    else:
        c = 2 * pyo.asin(pyo.sqrt(a + eps))
    # Radius of Earth in kilometers (use 3956 for miles)
    r = 6371
    return c * r

def linear_distance(lat1, lon1, lat2, lon2):
    lat_dist = lat_dist_conv(lat1, lat2)
    lon_dist = lon_dist_conv(lat1, lat2, lon1, lon2)
    return lat_dist + lon_dist

def euclidean_distance(lat1, lon1, lat2, lon2):
    lat_dist = lat_dist_conv(lat1, lat2)
    lon_dist = lon_dist_conv(lat1, lat2, lon1, lon2)
    eps = 1e-6  # to ensure that sqrt is never 0
    return pyo.sqrt(lat_dist**2+lon_dist**2+eps)

def lat_dist_conv(lat1, lat2):
    # from lat1, lat2 in degrees to a distance in km
    return 111.32 * abs(lat2-lat1)

def lon_dist_conv(lat1, lat2, lon1, lon2):
    # from lon1, lon2 in degrees to a distance in km (depends on lat1, lat2)
    avglat = (lat1+lat2)/2
    lon_dist = 111.32 * pyo.cos(avglat*math.pi/180) * abs(lon2 - lon1)
    return lon_dist