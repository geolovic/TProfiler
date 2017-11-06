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
#
#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com
#
#  Version: 2
#  29 October 2017

#  Last modified 30 October 2017

import ogr
import os
import numpy as np
import profiler as p
import pickle


# PROGRAM CODE
# =============
def save_profiles(path, profiles):
    """
    Saves the profiles in lists with pickle module
    Each profile will be a tuple with (_data, dem_res, rid, thetaref, chi0, slope_reg_points, _srs, name, _mouthdist, 0)

    :param path: Full path to the file where profiles will be saved
    :param profiles: List of TProfile objects
    """
    out_f = open(path, "wb")
    out_p = []
    for perfil in profiles:
        chi0 = perfil._data[perfil.n_points - 1, 5]
        out_p.append((perfil._data, perfil.dem_res, perfil.rid, perfil.thetaref, chi0, perfil.slope_reg_points,
                      perfil._srs, perfil.name, perfil._mouthdist, 0))

    pickle.dump(out_p, out_f)


def main(dem, fac, out_chi, out_file=False, threshold=0, units="CELL", basin_shp="", head_shp="",  id_field="",
         thetaref=0.45, reg_points=5, smooth=0, distance=0):

    # Get all the heads according to the given threshold
    if threshold:
        heads = p.get_heads(fac, dem, threshold, units)
    else:
        heads = np.array([], dtype="float32").reshape((0, 6))

    # Obtenemos todas las cabeceras de la capa de puntos, y combinamos con cabeceras anteriores
    if head_shp:
        main_heads = p.heads_from_points(dem, head_shp, id_field=id_field)
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
    if out_file:
        out_file = os.path.splitext(out_chi)[0] + ".dat"
        save_profiles(out_file, out_profiles)

    # Save output profiles in a shapefile
    p.profiles_to_shp(out_chi, out_profiles, distance)
