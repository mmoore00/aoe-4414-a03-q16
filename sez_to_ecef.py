# sez_to_ecef.py
#
# Usage: python3 sez_to_ecef.py lat_deg lon_deg hae_km s_km e_km z_km
# Converts SEZ coordinate vector to ECEF coordinate vector
# 
# Parameters:
#  lat_deg: latitude in degrees
#  lon_deg: longitude in degrees
#  hae_km: height above reference ellipsoid in kilometers
#  s_km: S distance value in SEZ frame in kilometers
#  e_km: E distance value in SEZ frame in kilometers
#  z_km: Z distance value in SEZ frame in kilometers
#  
# Output:
#  Print the ECEF coordinates
#
# Written by Matthew Moore
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# "constants"
R_E_KM = 6378.1363
E_E = 0.081819221456

# import Python modules
# e.g., import math # math module
import math # math module
import sys # argv

# helper functions

## matrix_times_vector
##
def matrix_times_vector(mat, vec, out_vec):
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            out_vec[i] += mat[i][j] * vec[j]
    return out_vec

## calc_denom
##
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-ecc**2.0 * math.sin(lat_rad)**2.0)
            

# initialize script arguments
rECEF_sez = [0, 0, 0]
Ry = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
Rz = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

# parse script arguments
if len(sys.argv)==7:
    lat_deg = float(sys.argv[1])
    lon_deg = float(sys.argv[2])
    hae_km = float(sys.argv[3])
    r_s_km = float(sys.argv[4])
    r_e_km = float(sys.argv[5])
    r_z_km = float(sys.argv[6])
else:
    print('Usage: python3 sez_to_ecef.py lat_deg lon_deg hae_km s_km e_km z_km ')
    exit()

# write script below this line

# initialize calculation variables
rSEZ = [r_s_km, r_e_km, r_z_km]
lat_rad = lat_deg * math.pi / 180
lon_rad = lon_deg * math.pi / 180

# initialize rotation matrices
Ry = [[math.sin(lat_rad), 0, math.cos(lat_rad)],
      [0, 1, 0],
      [-math.cos(lat_rad), 0, math.sin(lat_rad)]]
Rz = [[math.cos(lon_rad), -math.sin(lon_rad), 0],
      [math.sin(lon_rad), math.cos(lon_rad), 0],
      [0, 0, 1]]

# calculate Ry and Rz rotations
RyrSEZ = [0, 0, 0]
RyrSEZ = matrix_times_vector(Ry, rSEZ, RyrSEZ)
rECEF_sez = matrix_times_vector(Rz, RyrSEZ, rECEF_sez)

# calculate ECEF to SEZ origin
denom = calc_denom(E_E, lat_rad)
c_E = R_E_KM/denom
s_E = (R_E_KM * (1.0 - E_E * E_E)) / denom
r_x_km = (c_E + hae_km) * math.cos(lat_rad) * math.cos(lon_rad)
r_y_km = (c_E + hae_km) * math.cos(lat_rad) * math.sin(lon_rad)
r_z_km = (s_E + hae_km) * math.sin(lat_rad)

# find total ECEF vector
ecef_x_km = rECEF_sez[0] + r_x_km
ecef_y_km = rECEF_sez[1] + r_y_km
ecef_z_km = rECEF_sez[2] + r_z_km

print(rECEF_sez)
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)
