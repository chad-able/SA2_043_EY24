import modeling_functions
import cost_functions
import distance_functions

print(modeling_functions.sigmoid_shifted(0, 30, 1))

print(cost_functions.transportation_cost(0.1, 10))


h_approx = True

lat_0 = 59.12138156986758
lon_0 = -63.55275788848772
lat_1 = 5.506770208305376
lon_1 = -72.26282880715189

print("Haversine,", distance_functions.haversine(lat_0, lon_0, lat_1, lon_1, h_approx))
print("Linear,", distance_functions.linear_distance(lat_0, lon_0, lat_1, lon_1))
print("Euclidean,", distance_functions.euclidean_distance(lat_0, lon_0, lat_1, lon_1))