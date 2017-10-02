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

#  Last modified 02 October, 2017

import ogr
import os
import numpy as np
import argparse
import profiler as p

# ARGUMENT PARSER
# ===============
parser = argparse.ArgumentParser()
parser.add_argument("dem", help="Digital Elevation Model")
parser.add_argument("fac", help="Flow Accumulation Raster")
parser.add_argument("out_shp", help="Output Chi shapefile")
parser.add_argument("-th", "--threshold", help="Flow Accumulation threshold", type=float, default=0.)
parser.add_argument("-u", "--units",  help="Threshold units", choices=["CELL", "MAP"], default="CELL")
parser.add_argument("-b", "--basins", help="Basins shapefile", default="")
parser.add_argument("-hd", "--heads",  help="Heads shapefile", default="")
parser.add_argument("-id", "--id_field",  help="Id Field in heads shapefile", default="")
parser.add_argument("-t", "--thetaref", type=float, default=0.45, help="Reference n/m value for chi index evaluation")
dist_help = "Segment distance for output Chi shapefile. Distance = 0 will generate a point shapefile"
parser.add_argument("-d", "--distance", type=float, default=0., help=dist_help)
parser.add_argument("-r", "--reg_points", type=int, default=4, help="Number of points for slope and ksn regressions")
parser.add_argument("-s", "--smooth", type=float, default=0.0, help="Distance for smooth profile elevations")
args = parser.parse_args()


# ARGUMENTS
# ===============
dem = args.dem
fac = args.fac
threshold = args.threshold
units = args.units
basin_shp = args.basins
head_shp = args.heads
id_field = args.id_field
thetaref = args.thetaref
reg_points = args.reg_points
smooth = args.smooth
distance = args.distance
out_chi = args.out_shp
out_file = os.path.splitext(out_chi)[0] + ".npy"


# # DEBUG ARGUMENTS
# # ===============
# dem = "../../test/data/in/darro25.tif"
# fac = "../../test/data/in/darro25fac.tif"
# threshold = 1000
# units = "CELL"
# basin_shp = ""
# head_shp = ""
# id_field = ""
# thetaref = 0.45
# reg_points = 4
# smooth = 0
# distance = 200
# out_chi = "../../test/data/out/QGIS_ChiMap03.shp"
# out_file = "../../test/data/out/QGIS_ChiMap03.npy"


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

    # Save output profiles as a list of TProfiles
    np.save(out_file, np.array(out_profiles))

    # Save output profiles in a shapefile
    p.profiles_to_shp(out_chi, out_profiles, distance)


main(dem, fac, threshold, units, basin_shp, head_shp, id_field, thetaref, reg_points, smooth, distance, out_chi,
     out_file)
