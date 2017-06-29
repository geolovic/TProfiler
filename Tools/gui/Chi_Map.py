# -*- coding: iso-8859-15 -*-
#
#  Chi_Map.py
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

#  Version: 1.1
#  June 21, 2017

#  Last modified June 23, 2017

import ogr
import numpy as np
from profiler import get_heads, heads_from_points, heads_inside_basin, get_profiles, profiles_to_shp


# DEBUG ARGUMENTS
# ===============

basedir = "../../test/data/"
dem = basedir + "darro25.tif"
fac = basedir + "darro25fac.tif"
out_chi = basedir + "out_chi_map.shp"
out_file = basedir + "out_chi_profiles.npy"
distance = 250
threshold = 1000
units = "CELL"
basin_shp = basedir + "cuencas.shp"
head_shp = basedir + "main_heads.shp"
id_field = "id"
thetaref = 0.45
reg_points = 4
smooth = 250


# PROGRAM CODE
# =============


def main(dem, fac, out_chi, out_file, distance=0, umbral=0, units="CELL", basin_shp="", head_shp="",
         id_field=id_field, thetaref=0.45, reg_points=4, smooth=0):

    # Obtenemos todas las cabeceras del DEM
    if umbral:
        heads = get_heads(fac, dem, umbral, units)
    else:
        heads = np.array([], dtype="float32").reshape((0, 6))

    # Obtenemos los canales principales
    if head_shp:
        main_heads = heads_from_points(dem, head_shp, id_field)
    else:
        main_heads = np.array([], dtype="float32").reshape((0, 6))

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)
    if len(heads) == 0:
        return

    # Lista vacia para perfiles de salida
    out_profiles = []

    if basin_shp:
        dataset = ogr.Open(basin_shp)
        layer = dataset.GetLayer(0)
        for feat in layer:
            basin_geom = feat.GetGeometryRef()
            basin_heads = heads_inside_basin(heads, basin_geom)
            if len(basin_heads) > 0:
                out_profiles.extend(get_profiles(fac, dem, basin_heads, basin=basin_geom, tributaries=True,
                                                 thetaref=thetaref, reg_points=reg_points, smooth=smooth))
    else:
        out_profiles.extend(get_profiles(fac, dem, heads, tributaries=True, thetaref=thetaref, reg_points=reg_points,
                                         smooth=smooth))

    # Save output profiles
    np.save(out_file, np.array(out_profiles))

    # Store output profiles in a shapefile
    if distance > 0:
        profiles_to_shp(out_chi, out_profiles, distance)
    else:
        profiles_to_shp(out_chi, out_profiles)


main(dem, fac, out_chi, out_file, distance, threshold, units, basin_shp, head_shp, id_field, thetaref, reg_points,
     smooth)
