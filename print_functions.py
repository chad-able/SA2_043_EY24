import csv
import cost_functions
import folium
from folium.plugins import MarkerCluster

def print_initial(initial_cost, model, num_facilities, num_sites, site_coordinates, sites_flag):
    print_array = []

    print(f"Initial value for annualized cost: {initial_cost} $/bbl")
    col = []
    col.append(f"Initial value for annualized cost: {initial_cost} $/bbl")
    print_array.append(col)

    if not sites_flag:
        for i in range(num_facilities):
            col = []
            col.append(f"Facility {i}:")
            col.append(f"Location ({model.lat[i].value}, {model.lon[i].value})")
            print_array.append(col)
            print(f"Initial lat[{i}] = {model.lat[i].value}, lon[{i}] = {model.lon[i].value}")

    else:
        for i in range(num_facilities):
            lat, lon = cost_functions.facility_coordinates_vals(model, i, site_coordinates, num_sites)
            col = []
            col.append(f"Facility {i}:")
            col.append(f"Location ({lat}, {lon})")
            print_array.append(col)
            print(f"Initial lat[{i}] = {lat}, lon[{i}] = {lon}")
            for j in range(num_sites):
                if model.z[i, j].value == 1:
                    col = []
                    col.append(f"Facility {i} is initially assigned to Site {j}")
                    print_array.append(col)
                    print(f"Facility {i} is assigned to Site {j}")

    return print_array

def print_final(print_array, final_cost, model, num_facilities, num_sites, site_coordinates, flow_rate_data, sites_flag, filename):
    print_array.append("")
    print(f"Final value for annualized cost is {final_cost} $/bbl")
    col = []
    col.append(f"Final value for annualized cost is {final_cost} $/bbl")
    print_array.append(col)

    for i in range(num_facilities):
        if not sites_flag:
            col = []
            col.append(f"Facility {i}:")
            col.append(f"Location ({model.lat[i].value}, {model.lon[i].value})")
            print_array.append(col)
            print(f"Facility {i}: Location ({model.lat[i].value}, {model.lon[i].value})")
        else:
            lat, lon = cost_functions.facility_coordinates_vals(model, i, site_coordinates, num_sites)
            col = []
            col.append(f"Facility {i}:")
            col.append(f"Location ({lat}, {lon})")
            print_array.append(col)
            print(f"Facility {i}: Location ({lat}, {lon})")
        for j in range(num_sites):
            if model.x[i, j].value == 1:
                col = []
                col.append(f"Draws from site {j} with flow rate {flow_rate_data[j][0]}")
                print_array.append(col)
                print(f"  Draws from site {j} with flow rate {flow_rate_data[j][0]}")
            if sites_flag:
                if model.z[i, j].value == 1:
                    col = []
                    col.append(f"Facility {i} is assigned to Site {j}")
                    print_array.append(col)
                    print(f"Facility {i} is assigned to Site {j}")

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(print_array)):
            writer.writerow(print_array[i])



def mapping_flow_rates(shale_df, well_df, basin_stats):
    shale_df = shale_df.merge(basin_stats, on="Shale_play", how="left")
    well_df = well_df.drop(columns=["loaddate", "geometry_centroids"])
    basin_wgs84 = shale_df.to_crs(epsg=4326)
    m = folium.Map(location=[39.5, -98.35], zoom_start=5)
    folium.GeoJson(basin_wgs84, tooltip=folium.features.GeoJsonTooltip(fields=['Shale_play'], aliases=['Basin:'])).add_to(m)
    marker_group = folium.FeatureGroup(name="HUC8 Markers", show=True)
    for idx, row in well_df.iterrows():

        lat, lon = row["lat"], row["lon"]
        flow_rate = row["2022_flow_gpm"] if "2022_flow_gpm" in row else "N/A"

        huc8_id = row["huc8"] if "huc8" in row else "N/A"

        popup_text = f"HUC-8: {huc8_id}<br>Flow Rate: {flow_rate}<br>Coords: ({lat:.3f}, {lon:.3f})"
        folium.Marker(location=[lat, lon], popup=popup_text).add_to(marker_group)

    marker_group.add_to(m)
    folium.LayerControl().add_to(m)
    folium.GeoJson(well_df).add_to(m)
    m.save("map.html")
