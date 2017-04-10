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
    Cuenca del Darro
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
    out_txt = "data/02_basins_test01" + ".txt"

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
    """
    Test for profiler.heads_inside_basin() function
    Cuencas.shp
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print "=" * 40
    print "Test 02 para profiler.heads_inside_basin() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_channels.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_test02" + ".txt"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos el poligono de la tercera feature (Aguas blancas)
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(2)
    geom = feat.GetGeometryRef()

    # Obtenemos las cabeceras de dentro de la cuenca
    basin_heads = p.heads_inside_basin(heads, geom)

    # Obtenemos los canales principales dentro del poligono
    main_heads = p.heads_from_points(dem, main_ch)
    main_heads = p.heads_inside_basin(main_heads, geom)

    # Combinamos ambos arrays (anadimos las cab principales al principio)
    for head in main_heads:
        basin_heads.insert(0, head)

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


def test03():
    """
    Test for profiler.heads_inside_basin() function
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print "=" * 40
    print "Test 03 para profiler.heads_inside_basin() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_channels.shp"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_test03" + ".txt"

    # Creamos archivo de salida y escribimos primera linea
    outfile = open(out_txt, "w")
    outfile.write("id;ROW;COL;X;Y;Z;Name\n")

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)
    # Obtenemos los canales principales
    main_heads = p.heads_from_points(dem, main_ch)

    # Obtenemos los diferentes poligonos de las cuencas
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)

    id_basin = 0

    for feat in layer:
        geom = feat.GetGeometryRef()
        # Obtenemos las cabeceras de dentro de la cuenca
        basin_heads = p.heads_inside_basin(heads, geom)
        basin_main_heads = p.heads_inside_basin(main_heads, geom)

        # Combinamos ambos arrays (anadimos las cab principales al principio)
        for head in basin_main_heads:
            basin_heads.insert(0, head)

        # Escribimos en archivo de texto
        for head in basin_heads:
            head = [str(value) for value in head]
            linea = ";".join(head) + "\n"
            outfile.write(str(id_basin) + ";" + linea)

        id_basin += 1
    outfile.close()

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "Resultado en " + out_txt
    print "=" * 40


test01()
test02()
test03()
