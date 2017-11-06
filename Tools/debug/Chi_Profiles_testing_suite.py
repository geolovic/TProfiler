# -*- coding: iso-8859-15 -*-
#
#  Testing Suite for >> Chi_Profiles.py
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

from Chi_Profiles import main
import time


def test_00():
    """
    Test 00 for Chi_Profiles
    Ejecuta Chi_Profiles con un shapefile de rios

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 00 para Chi_Profiles")
    print("Executes Chi_Profiles for a river shapefile")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    river_shp = "../data/in/rios.shp"
    out_file = "../data/out/perfiles_rios.npy"
    id_field = "id"
    name_field = "name"
    thetaref = 0.45
    reg_points = 10
    smooth = 0

    main(dem, fac, river_shp, out_file, id_field=id_field, name_field=name_field, thetaref=thetaref,
         reg_points=reg_points, smooth=smooth)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output in " + out_file)
    print("=" * 40)


def test_01():
    """
    Test 01 for Chi_Profiles
    Ejecuta Chi_Profiles con un shapefile de rios que no están conectados

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 00 para Chi_Profiles")
    print("Executes Chi_Profiles for a river shapefile")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    river_shp = "../data/in/rios2.shp"
    out_file = "../data/out/Chi_Profiles_test01.npy"
    id_field = "id"
    name_field = "name"
    thetaref = 0.45
    reg_points = 10
    smooth = 0

    main(dem, fac, river_shp, out_file, id_field=id_field, name_field=name_field, thetaref=thetaref,
         reg_points=reg_points, smooth=smooth)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output in " + out_file)
    print("=" * 40)


test_00()
test_01()
