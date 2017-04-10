# -*- coding: iso-8859-15 -*-
#
#  Chi_Map_import.py
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

#  Version: 1.0
#  April 08, 2017

#  Last modified April 08, 2017

import ogr
import profiler as p

# DEM=raster
# Flow_Accumulation=raster
# Output_shapefile=output vector
# Threshold=number 1000
# Units=selection CELLS;MAP
# Distance=number 0
# Basins=boolean False
# Basin_shapefile=vector
# Main_channels=boolean False
# Head_shapefile=vector
# Thetaref=number 0.45
# Regression_points=number 4
# Smooth=number 0

# Debug code
DEM = "../../test/data/darro25.tif"
Flow_Accumulation = "../../test/data/darro25fac.tif"
Distance = 250
Output_shapefile = "../../test/data/qgis_test_import.shp"
Threshold = 1000
Units = "CELL"
Basins = True
Basin_shapefile = "../../test/data/cuencas.shp"
Main_channels = True
Head_shapefile = "../../test/data/main_channels.shp"
Thetaref = 0.45
Regression_points = 4
Smooth = 0


def main(dem, fac, distance, out_shapefile, umbral, units, basins, basin_shp, main_ch, head_shp, thetaref,
         reg_points, smooth):

    if not basins:
        basin_shp = ""
    if not main_ch:
        head_shp = ""

    # Get all the heads of the dem and incorporate main channels if specified
    heads = p.get_heads(fac, dem, umbral, units)
    if main_ch:
        main_heads = p.heads_from_points(dem, head_shp)
        main_heads.reverse()
        for head in main_heads:
            heads.insert(0, head)

    # Empty list to store output profiles
    out_profiles = []

    if basin_shp:
        dataset = ogr.Open(basin_shp)
        layer = dataset.GetLayer(0)
        for feat in layer:
            basin_geom = feat.GetGeometryRef()
            basin_heads = p.heads_inside_basin(heads, basin_geom)
            if len(basin_heads) > 0:
                out_profiles.extend(p.get_profiles(fac, dem, basin_heads, basin=basin_geom, tributaries=True,
                                                   thetaref=thetaref, reg_points=reg_points, smooth=smooth))
    else:
        out_profiles.extend(p.get_profiles(fac, dem, heads, tributaries=True))

    # Store output profiles in a shapefile
    if distance > 0:
        p.profiles_to_shp(out_shapefile, out_profiles, distance)
    else:
        p.profiles_to_shp(out_shapefile, out_profiles)


main(DEM, Flow_Accumulation, Distance, Output_shapefile, Threshold, Units, Basins,
     Basin_shapefile, Main_channels, Head_shapefile, Thetaref, Regression_points, Smooth)
