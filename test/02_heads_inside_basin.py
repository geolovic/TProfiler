# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 09 March 2017
"""

import time
import profiler as p

print "Test for profiler.heads_inside_basin()"

def test01():
    """
    Test for profiler.heads_inside_basin() function
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print "=" * 40
    print "Test 01 para profiler.heads_inside_basin() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuenca_darro.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_darro_heads" + str(umbral) + units + ".txt"

    cabeceras = p.heads_inside_basin(fac, dem, basin, umbral, units)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "Resultado en " + out_txt
    print "=" * 40


def test02():
    """
    Test for profiler.heads_inside_basin() function
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print "=" * 40
    print "Test 02 para profiler.heads_inside_basin() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuenca_darro.shp"
    main_ch = "data/darro_main.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_darro_main_head" + str(umbral) + units + ".txt"

    cabeceras = p.heads_inside_basin(fac, dem, basin, umbral, units, main_ch)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

    for cab in cabeceras:
        cab = [str(value) for value in cab]
        linea = ";".join(cab) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "Resultado en " + out_txt
    print "=" * 40


test01()
test02()
