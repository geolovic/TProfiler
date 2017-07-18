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
import profiler as p


# ARGUMENTS
# ===============
dem = "../../test/data/in/darro25.tif"
fac = "../../test/data/in/darro25fac.tif"
threshold = 1000
units = "CELL"
basin_shp = "../../test/data/in/cuencas.shp"
head_shp = "../../test/data/in/main_heads.shp"
id_field = "id"
thetaref = 0.45
reg_points = 4
smooth = 0
distance = 0
out_chi = "../../test/data/out/ChiMap_output_shp.shp"
out_file = "../../test/data/out/ChiMap_output_file.npy"


# PROGRAM CODE
# =============
def main(dem, fac, threshold, units, basin_shp, head_shp, id_field, thetaref, reg_points, smooth, distance, out_chi,
         out_file):

    # Obtenemos todas las cabeceras del DEM
    if threshold:
        heads = p.get_heads(fac, dem, threshold, units)
    else:
        heads = np.array([], dtype="float32").reshape((0, 6))

    # Obtenemos los canales principales
    if head_shp:
        main_heads = p.heads_from_points(dem, head_shp, id_field)
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
            heads_inside = p.heads_inside_basin(heads, basin_geom)
            if len(heads_inside) > 0:
                out_profiles.extend(p.get_profiles(fac, dem, heads_inside, basin=basin_geom, tributaries=True,
                                                 thetaref=thetaref, reg_points=reg_points, smooth=smooth))
    else:
        out_profiles.extend(p.get_profiles(fac, dem, heads, tributaries=True, thetaref=thetaref, reg_points=reg_points,
                                         smooth=smooth))

    # Save output profiles
    np.save(out_file, np.array(out_profiles))

    # Save output profiles in a shapefile
    p.profiles_to_shp(out_chi, out_profiles, distance)


main(dem, fac, threshold, units, basin_shp, head_shp, id_field, thetaref, reg_points, smooth, distance, out_chi,
     out_file)
