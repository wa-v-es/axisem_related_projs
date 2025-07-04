import math
# Constants
earth_radius_km = 6371 # Earth's radius in kilometers
station_spacing_km = 200 # Approximate spacing between stations on each ring
start_radius_deg = 4
end_radius_deg = 45
# Function to convert degrees to radians
def deg2rad(deg):
    return deg * math.pi / 180
# Function to generate station names
def generate_station_name(index):
    letters = ''
    while index >= 0:
        letters = chr(index % 26 + ord('A')) + letters
        index = index // 26 - 1
    return letters

# Generate stations
stations = []
station_index = 0

for radius_deg in range(start_radius_deg, end_radius_deg + 1):
    radius_rad = deg2rad(radius_deg)
    arc_length_km = radius_rad * earth_radius_km
    num_stations = max(1, round(2 * math.pi * arc_length_km / station_spacing_km))

    for i in range(num_stations):
        angle_rad = 2 * math.pi * i / num_stations
        lat = radius_deg * math.cos(angle_rad)
        lon = radius_deg * math.sin(angle_rad)
        name_prefix = generate_station_name(station_index)
        station_name = f"{name_prefix}{i+1}"
        stations.append((station_name, 'II', lat, lon, 0.0, 0.0))

    station_index += 1

with open("stations_output.txt", "w") as f:
    f.write("# ------------------------------------------------------------------------------\n")
    f.write("#     name   network            latitude           longitude  useless      depth\n")
    f.write("# ------------------------------------------------------------------------------\n")
    for station in stations:
        f.write(f"{station[0]:>8}   {station[1]:>4}   {station[2]:>15.6f}   {station[3]:>15.6f}   {station[4]:>7.1f}   {station[5]:>7.1f}\n")

print("Station data has been written to 'stations_output.txt'.")
