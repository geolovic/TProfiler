# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 16 June 2017
"""

import time
import profiler as p

print("Test for profiler.heads_from_points()")


def test01():
    """
    Test for profiler.heads_from_points() function
    Intput points are in "data/in/cabeceras_darro.shp"
    """
    
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiler.heads_from_points() function")
    print("Testing a right id_field")
    print("Test in progress...")
    
    # Test parameters
    dem = "data/in/darro25.tif"
    pointshp = "data/in/cabeceras_darro.shp"
    out_txt = "data/out/01_cabeceras_darro_1.txt"
    id_field = "id"
    
    cabeceras = p.heads_from_points(dem, pointshp, id_field)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;id\n")
    
    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)
    
    outfile.close()
    
    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


def test02():
    """
    Test for profiler.heads_from_points() function
    Intput points are in "data/in/cabeceras_darro.shp"
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 02 para profiler.heads_from_points() function")
    print("Testing without id_field")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    pointshp = "data/in/cabeceras_darro.shp"
    out_txt = "data/out/01_cabeceras_darro_2.txt"

    cabeceras = p.heads_from_points(dem, pointshp)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;id\n")

    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


def test03():
    """
    Test for profiler.heads_from_points() function
    Intput points are in "data/in/cabeceras_darro.shp"
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 03 para profiler.heads_from_points() function")
    print("Testing a field that is not in head shapefile")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    pointshp = "data/in/cabeceras_darro.shp"
    out_txt = "data/out/01_cabeceras_darro_3.txt"
    id_field = "id_rio"

    cabeceras = p.heads_from_points(dem, pointshp, id_field)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;id\n")

    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


def test04():
    """
    Test for profiler.heads_from_points() function
    Intput points are in "data/in/cabeceras_darro.shp"
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 04 para profiler.heads_from_points() function")
    print("Testing wrong id_field (bad field type)")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    pointshp = "data/in/cabeceras_darro.shp"
    out_txt = "data/out/01_cabeceras_darro_4.txt"
    id_field = "Name"

    cabeceras = p.heads_from_points(dem, pointshp, id_field)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;id\n")

    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)

test01()
test02()
test03()
test04()
