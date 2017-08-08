# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 16 June 2017
    - Adapted for Python 3.5
"""

import time
import profiler as p
print("-" * 40)
print("Tests for profiler.get_heads()")


def test01():
    """
    Test for get_heads() function
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiler.get_heads() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    umbral = 1000
    units = "CELL"
    out_txt = "data/out/00_heads_" + str(umbral) + units + ".txt"

    cabeceras = p.get_heads(fac, dem, umbral, units)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

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
    Test for get_heads() function
    Umbral de 625000 m2 (1000 celdas)
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 02 para profiler.get_heads() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    umbral = 625000
    units = "MAP"
    out_txt = "data/out/00_heads_" + str(umbral) + units + ".txt"

    cabeceras = p.get_heads(fac, dem, umbral, units)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

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
    Test for get_heads() function
    Umbral de 1000 celdas (sin especificar units)
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 03 para profiler.get_heads() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    umbral = 1000
    out_txt = "data/out/00_heads_" + str(umbral) + ".txt"

    cabeceras = p.get_heads(fac, dem, umbral)
    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

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
