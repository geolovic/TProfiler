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

#  Version: 1.1
#  July 18, 2017

#  Last modified 29 October, 2017

import numpy as np
import profiler as p


# PROGRAM CODE
# =============
def main(dem, fac, river_shp, out_file, id_field="", name_field="", thetaref=0.45, reg_points=4, smooth=0):
    # Extract profiles
    perfiles = p.profiles_from_rivers(fac, dem, river_shp, id_field=id_field, name_field=name_field, thetaref=thetaref,
                                      reg_points=reg_points, smooth=smooth)
    # Save profiles into file
    perfiles = np.array(perfiles)
    np.save(out_file, perfiles)
