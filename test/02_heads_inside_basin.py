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
import numpy as np
import ogr

print("Test for profiler.heads_inside_basin()")


def test01():
    """
    Test for profiler.heads_inside_basin() function
    Cuenca del Darro
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiler.heads_inside_basin() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuenca_darro.shp"
    cabeceras = "data/cabeceras_darro.shp"
    id_field = "id"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_test01" + ".txt"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras del Darro
    main_heads = p.heads_from_points(dem, cabeceras, id_field=id_field)

    # Combinamos las cabeceras del Darro con las demas del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos el poligono de la primera feature (cuenca Darro)
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(0)
    geom = feat.GetGeometryRef()

    # Obtenemos las cabeceras dentro de la cuenca
    basin_heads = p.heads_inside_basin(heads, geom)
    basin_heads[:, 5] = np.arange(basin_heads.shape[0])

    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Id\n")

    for head in basin_heads:
        head = [str(value) for value in head]
        linea = ";".join(head) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


def test02():
    """
    Test for profiler.heads_inside_basin() function
    Cuencas.shp
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 02 para profiler.heads_inside_basin() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_heads.shp"
    id_field = "id"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_test02" + ".txt"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_ch, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos el poligono de la tercera feature (Aguas blancas)
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(2)
    geom = feat.GetGeometryRef()

    # Obtenemos las cabeceras dentro de la cuenca del Aguas Blancas
    basin_heads = p.heads_inside_basin(heads, geom)
    basin_heads[:, 5] = np.arange(basin_heads.shape[0])
    print(basin_heads.dtype)

    outfile = open(out_txt, "w")
    outfile.write("ROW;COL;X;Y;Z;Id\n")

    for head in basin_heads:
        head = [str(value) for value in head]
        linea = ";".join(head) + "\n"
        outfile.write(linea)

    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


def test03():
    """
    Test for profiler.heads_inside_basin() function
    Umbral de 1000 celdas
    """

    inicio = time.time()
    print("=" * 40)
    print("Test 03 para profiler.heads_inside_basin() function")
    print("Test in progress...")

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_heads.shp"
    id_field = "id"
    umbral = 1000
    units = "CELL"
    out_txt = "data/02_basins_test03" + ".txt"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_ch, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos los diferentes poligonos de las cuencas
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)

    id_basin = 0

    # Creamos archivo de salida y escribimos primera linea
    outfile = open(out_txt, "w")
    outfile.write("Basin;ROW;COL;X;Y;Z;Id\n")

    n_heads = 0

    for feat in layer:

        # Obtenemos la cuenca y las cabeceras dentro de la misma
        geom = feat.GetGeometryRef()
        basin_heads = p.heads_inside_basin(heads, geom, n_heads)
        n_heads += basin_heads.shape[0]

        # Escribimos en archivo de texto
        for head in basin_heads:
            head = [str(value) for value in head]
            linea = ";".join(head) + "\n"
            outfile.write(str(id_basin) + ";" + linea)

        id_basin += 1
    outfile.close()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Resultado en " + out_txt)
    print("=" * 40)


#test01()
test02()
#test03()
