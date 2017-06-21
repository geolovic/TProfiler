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

#  Last modified June 21, 2017

import ogr
import numpy as np

from profiler import get_heads, heads_from_points, heads_inside_basin, get_profiles, profiles_to_shp, TProfile

##DEM=raster
##Flow_accumulation=raster
##Output_chi=output vector
##Use_threshold=boolean True
##Threshold=number 1000
##Units=selection CELLS;MAP
##Use_basins=boolean True
##Basin_shapefile=vector
##Use_heads=boolean True
##Heads_shapefile=vector
##Id_field=field Heads_shapefile
##Thetaref=number 0.45
##Regression_points=number 4
##Smooth=number 0

# Debug code
basedir = "../../test/data/"
DEM = basedir + "darro25.tif"
Flow_Accumulation = basedir + "darro25fac.tif"
Use_threshold = True
Threshold = 1000
Units = "CELL"
Use_basins = True
Basin_shapefile = basedir + "cuencas.shp"
Use_heads = True
Head_shapefile = basedir + "main_heads.shp"
Id_field = "id"
Thetaref = 0.45
Regression_points = 4
Smooth = 0
Distance = 250
Output_chi = basedir + "out_chi_map.shp"
# End Debug


def main(dem, fac, distance, out_shapefile, use_umbral, umbral, units, use_basins, basins_shp, use_heads, heads_shp,
         thetaref, reg_points, smooth):

    # Obtenemos todas las cabeceras del DEM
    if use_umbral:
        heads = get_heads(fac, dem, umbral, units)
    else:
        heads = np.array([], dtype="float32").reshape((0, 6))

    # Obtenemos los canales principales
    if use_heads:
        main_heads = heads_from_points(dem, heads_shp)
    else:
        main_heads = np.array([], dtype="float32").reshape((0, 6))

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Lista vacia para perfiles de salida
    out_profiles = []

    if use_basins:
        dataset = ogr.Open(basins_shp)
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

    # Store output profiles in a shapefile
    if distance > 0:
        profiles_to_shp(out_shapefile, out_profiles, distance)
    else:
        profiles_to_shp(out_shapefile, out_profiles)


main(DEM, Flow_Accumulation, Distance, Output_chi, Use_threshold, Threshold, Units, Use_basins,
     Basin_shapefile, Use_heads, Head_shapefile, Thetaref, Regression_points, Smooth)
