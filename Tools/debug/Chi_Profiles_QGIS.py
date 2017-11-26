# -*- coding: iso-8859-15 -*-
#
#  Chi_Profiles.py [QGIS Version]
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
#  November 06, 2017

#  Last modified 06 November, 2017

import profiler as p
import pickle
import os


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


def main(dem, fac, river_shp, out_file, id_field="", name_field="", thetaref=0.45, reg_points=4, smooth=0):
    # Extract profiles
    perfiles = p.profiles_from_rivers(fac, dem, river_shp, id_field=id_field, name_field=name_field, thetaref=thetaref,
                                      reg_points=reg_points, smooth=smooth)

    # Extension for the output file has to be ".dat"
    if not os.path.splitext(out_file)[1] == ".dat":
        out_file = os.path.splitext(out_file)[0] + ".dat"

    # Save profiles into file
    save_profiles(out_file, perfiles)
