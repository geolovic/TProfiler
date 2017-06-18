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
    Intput points are in "data/in/cabeceras.shp"
    """
    
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiler.heads_from_points() function")
    print("Test in progress...")
    
    # Test parameters
    dem = "data/darro25.tif"
    pointshp = "data/cabeceras.shp"
    out_txt = "data/01_cabeceras.txt"
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


test01()
