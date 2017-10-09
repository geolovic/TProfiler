# -*- coding: iso-8859-15 -*-
#
#  Profiles_to_shapefile.py
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

#  Version: 1
#  July 18, 2017

#  Last modified July 18, 2017

import numpy as np
from profiler import profiles_to_shp


# ARGUMENTS
# ===============
in_profiles = "../../test/data/in/profiles_basins.npy"
out_shapefile = "../../test/data/out/Profiles_To_Shapefile_test.shp"
distance = 0
reg_points = 4


# PROGRAM CODE
# =============
def main(in_profiles, out_shapefile, distance=0, reg_points=4):
    # Load profiles
    perfiles = np.load(in_profiles)

    # Change reg_points if output shapefile is a point shapefile
    if distance == 0:
        for perfil in perfiles:
            perfil.calculate_slope(reg_points=reg_points)
            perfil.calculate_ksn(reg_points=reg_points)

    # Store profiles in a shapefile
    profiles_to_shp(out_shapefile, perfiles, distance)


main(in_profiles, out_shapefile, distance, reg_points)
