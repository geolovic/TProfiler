# -*- coding: iso-8859-15 -*-
#
#  Chi_Profiles.py
#
#  Copyright (C) 2017  J. Vicente Perez, Universidad de Granada
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com

#  Version: 1.2
#  November, 6th 2017

#  Last modified November, 6th 2017


import numpy as np
import profiler as p
import argparse


# ARGUMENT PARSER
# ===============
parser = argparse.ArgumentParser()
parser.add_argument("dem", help="Digital Elevation Model")
parser.add_argument("fac", help="Flow Accumulation Raster")
parser.add_argument("rivers", help="River shapefile")
parser.add_argument("out_profiles", help="Output profiles npy file")
parser.add_argument("-id", "--id_field",  help="Field with id numbers in river shapefile", default="")
parser.add_argument("-n", "--name_field",  help="Field with names in river shapefile", default="")
parser.add_argument("-t", "--thetaref", type=float, default=0.45, help="Reference n/m value for chi index evaluation")
parser.add_argument("-r", "--reg_points", type=int, default=4, help="Number of points for slope and ksn regressions")
parser.add_argument("-s", "--smooth", type=float, default=0.0, help="Distance for smooth profile elevations")
args = parser.parse_args()


# ARGUMENTS
# ===============
dem = args.dem
fac = args.fac
river_shp = args.rivers
id_field = args.id_field
name_field = args.name_field
out_file = args.out_profiles
thetaref = args.thetaref
reg_points = args.reg_points
smooth = args.smooth


# PROGRAM CODE
# =============
def main(dem, fac, river_shp, out_file, id_field="", name_field="", thetaref=0.45, reg_points=4, smooth=0):
    # Extract profiles
    perfiles = p.profiles_from_rivers(fac, dem, river_shp, id_field=id_field, name_field=name_field, thetaref=thetaref,
                                      reg_points=reg_points, smooth=smooth)
    # Save profiles into file
    perfiles = np.array(perfiles)
    np.save(out_file, perfiles)

main(dem, fac, river_shp, out_file, id_field, name_field, thetaref, reg_points, smooth)
