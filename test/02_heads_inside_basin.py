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
import ogr

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
    main_ch = "data/darro_main.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_darro_heads" + str(umbral) + units + ".txt"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos la cabecera del canal principal
    main_head = p.heads_from_points(dem, main_ch)

    # Obtenemos el poligono de la primera feature (cuenca Darro)
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(0)
    geom = feat.GetGeometryRef()

    # Obtenemos las cabeceras de dentro
    basin_heads = p.heads_inside_basin(heads, geom)

    # Anadimos la cabecera del canal principal
    basin_heads.insert(0, main_head[0])

    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Name\n")

    for head in basin_heads:
        head = [str(value) for value in head]
        linea = ";".join(head) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "Resultado en " + out_txt
    print "=" * 40


def test02():
    # todo Revisar test, funciona mal
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
    cuencas = "data/cuencas.shp"
    main_chs = "data/main_channels.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_heads" + str(umbral) + units + ".txt"

    # Archivo para escribir output
    outfile = open(out_txt, "w")
    outfile.write("cuenca;ROW;COL;X;Y;Z;Name\n")

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos las cabeceras de los canales principales
    main_heads = p.heads_from_points(dem, main_chs)

    # Obtenemos los polígonos de las cuencas y sus cabeceras recursivamente
    dataset = ogr.Open(cuencas)
    layer = dataset.GetLayer(0)
    id_cuenca = 0
    for feat in layer:
        feat = layer.GetFeature(0)
        geom = feat.GetGeometryRef()
        cabeceras = p.heads_inside_basin(heads, geom)
        cab_basin = p.heads_inside_basin(main_heads, geom)
        for head in cab_basin:
            cabeceras.insert(0, head)

        for cab in cabeceras:
            cab = [str(value) for value in cab]
            linea = str(id_cuenca) + ";" + ";".join(cab) + "\n"
            outfile.write(linea)
        id_cuenca += 1

    outfile.close()

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "Resultado en " + out_txt
    print "=" * 40


test01()
test02()
