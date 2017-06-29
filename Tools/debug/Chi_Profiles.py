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

#  Version: 1.0
#  June 23, 2017

#  Last modified June 23, 2017

import numpy as np
from profiler import profiles_from_rivers


# DEBUG ARGUMENTS
# ===============

basedir = "../../test/data/"
dem = basedir + "darro25.tif"
fac = basedir + "darro25fac.tif"
river_shp = basedir + "rios.shp"
id_field = "id"
name_field = "name"
out_file = basedir + "river_chi_profiles.npy"
thetaref = 0.45
reg_points = 4
smooth = 250

# PROGRAM CODE
# =============


NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}


def main(dem, fac, river_shp, out_file, id_field="", name_field="", thetaref=thetaref, reg_points=reg_points,
         smooth=smooth):
    # Extract profiles
    perfiles = profiles_from_rivers(fac, dem, river_shp, id_field=id_field, name_field=name_field, thetaref=thetaref,
                                    reg_points=reg_points, smooth=smooth)
    # Save profiles into
    perfiles = np.array(perfiles)
    np.save(out_file, perfiles)


main(dem, fac, river_shp, out_file, id_field, name_field, thetaref, reg_points, smooth)
