# -*- coding: iso-8859-15 -*-
#
#  Chi_Profiles_import.py
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
#  June 20, 2017

#  Last modified June 20, 2017

import numpy as np
from profiler import profiles_from_rivers

##DEM=raster
##Flow_accumulation=raster
##River_shapefile=vector
##Id_field=field River_shapefile
##Name_field=field River_shapefile
##Output_profiles=output file


NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}


def main(dem, fac, river_shapefile, id_field, name_field, out_profiles):
    # Extract profiles
    perfiles = profiles_from_rivers(fac, dem, river_shapefile, id_field=id_field, name_field=name_field)
    # Save profiles into a txt file
    perfiles = np.array(perfiles)
    np.save(out_profiles, perfiles)

# Imported modules and functions [version June 20, 2017]



# Debug
DEM = "../../test/data/darro25.tif"
Flow_accumulation = "../../test/data/darro25fac.tif"
River_shapefile = "../../test/data/rios.shp"
Id_field = "id"
Name_field = "name"
Output_profiles = "../../test/data/chi_profiles.npy"
# End debug


main(DEM, Flow_accumulation, River_shapefile, Id_field, Name_field, Output_profiles)
