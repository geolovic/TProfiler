# -*- coding: iso-8859-15 -*-
#
#  Testing Suite for >> Get_Channels.py
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

from Get_Channels import main
import time


def test_00():
    """
    Test 00 for Get_Channels
    Ejecuta Get_Channels sin ningun argumento opcional

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 00 para Get_Channels()")
    print("Executes Get_Channels without optional arguments")
    print("No threshold, No heads, No basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test00.shp"

    main(dem, fac, out_shp)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("No se debe generar shapefile al no haber cabeceras")
    print("=" * 40)


def test_01():
    """
    Test 01 for Get_Channels
    Ejecuta Get_Channels utilizando solamente un umbral en Celdas

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para Get_Channels()")
    print("Executes Get_Channels for a threshold of 1000 CELL")
    print("No heads, No basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test01.shp"
    threshold = 1000
    units = "CELL"

    main(dem, fac, out_shp, threshold=threshold, units=units)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_shp)
    print("=" * 40)


def test_02():
    """
    Test 02 for Get_Channels
    Ejecuta Get_Channels utilizando solamente un umbral en Map Units

    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para Get_Channels()")
    print("Executes Get_Channels for a threshold of 625000 m^2")
    print("No heads, No basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test02.shp"
    threshold = 625000
    units = "MAP"

    main(dem, fac, out_shp, threshold=threshold, units=units)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_shp)
    print("=" * 40)


def test_03():
    """
    Test 03 for Get_Channels
    Ejecuta Get_Channels utilizando solo las cabeceras principales
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 03 para Get_Channels()")
    print("Executes Get_Channels only for selected heads")
    print("No Threshold, No basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test03.shp"
    head_shp = "../data/in/main_heads.shp"
    id_field = "river_id"

    main(dem, fac, out_shp, head_shp=head_shp, id_field=id_field)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_shp)
    print("=" * 40)


def test_04():
    """
    Test 04 for Get_Channels
    Ejecuta Get_Channels utilizando cabeceras y un umbral
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 04 para Get_Channels()")
    print("Executes Get_Channels for selected heads and a threshold")
    print("No basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test04.shp"
    head_shp = "../data/in/main_heads.shp"
    id_field = "river_id"

    main(dem, fac, out_shp, threshold=1000, units="CELL", head_shp=head_shp, id_field=id_field)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_shp)
    print("=" * 40)


def test_05():
    """
    Test 05 for Get_Channels
    Ejecuta Get_Channels utilizando cabeceras, un umbral y una capa de cuencas
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 05 para Get_Channels()")
    print("Executes Get_Channels for selected heads, using a threshold, and inside basins")
    print("Test in progress...")

    # Test parameters
    # ===============
    dem = "../data/in/darro25.tif"
    fac = "../data/in/darro25fac.tif"
    out_shp = "../data/out/GetChannels_test05.shp"
    head_shp = "../data/in/main_heads.shp"
    id_field = "river_id"
    basin_shp = "../data/in/cuencas.shp"

    main(dem, fac, out_shp, threshold=1000, units="CELL", basin_shp=basin_shp, head_shp=head_shp, id_field=id_field)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_shp)
    print("=" * 40)


test_00()
test_01()
test_02()
test_03()
test_04()
test_05()
