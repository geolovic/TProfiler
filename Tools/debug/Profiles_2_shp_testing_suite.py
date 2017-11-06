# -*- coding: iso-8859-15 -*-
#
#  Testing Suite for >> Profiles_2_shp.py
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

from Profiles_2_shp import main
import time


def test_01():
    """
    Test 00 for Profiles_2_shp
    Output normal

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para Profiles_2_shp()")
    print("Save as point shapefile and changing reg_points")
    print("Input file does not come from QGIS")
    print("Test in progress...")

    # Test parameters
    # ===============
    in_profiles = "../data/in/tprofiles_t04.npy"
    out_shapefile = "../data/out/Profiles2shp_test01.shp"
    reg_points = 10

    main(in_profiles, out_shapefile, reg_points=reg_points)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Otput file en " + out_shapefile)
    print("=" * 40)


def test_02():
    """
    Test 02 for Get_Channels
    Output de QGIS

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para Profiles_2_shp")
    print("Save as point shapefile and changing reg_points and distance")
    print("Input file comes from QGIS")
    print("Test in progress...")

    # Test parameters
    # ===============
    in_profiles = "../data/in/tprofiles_t04.dat"
    out_shapefile = "../data/out/Profiles2shp_test02.shp"
    reg_points = 10
    distance = 200

    main(in_profiles, out_shapefile, qgis=True, reg_points=reg_points, distance=distance)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Otput file en " + out_shapefile)
    print("=" * 40)


test_01()
test_02()
