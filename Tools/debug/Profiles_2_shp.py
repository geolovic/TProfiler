# -*- coding: iso-8859-15 -*-
#
#  Profiles_2_shp.py
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

#  Last modified 30 October, 2017

import numpy as np
import pickle
from profiler import profiles_to_shp, TProfile


# PROGRAM CODE
# =============
def load_qgis_data(in_file):
    profile_data = pickle.load(open(in_file, "rb"), encoding="latin1")
    perfiles = []
    for row in profile_data:
        perfiles.append(TProfile(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], 0))
    return perfiles


def main(in_profiles, out_shapefile, distance=0, qgis=False, reg_points=0):
    # Load profiles
    if qgis:
        perfiles = load_qgis_data(in_profiles)
    else:
        perfiles = np.load(in_profiles)

    # Change reg_points if output shapefile is a point shapefile
    if distance == 0 and reg_points > 0:
        for perfil in perfiles:
            perfil.calculate_slope(reg_points=reg_points)
            perfil.calculate_ksn(reg_points=reg_points)

    # Store profiles in a shapefile
    profiles_to_shp(out_shapefile, perfiles, distance)
