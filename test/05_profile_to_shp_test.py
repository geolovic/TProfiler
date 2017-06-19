# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
March, 2017
Testing suite for profiler.py
Last modified: 19 June 2017
"""

import time
import profiler as p
import ogr
import numpy as np

print("Tests for profiler.profiles_to_shp()")

    
def test01():
    """
    Test for profiles_to_shp() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_to_shp() function")
    print("Output point and line shapefiles")
    print("Test in progress...")
    
    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_heads.shp"
    umbral = 1000
    units = "CELL"
    out_point_shp = "data/05_profile_points" + ".shp"
    dl = 250
    out_line_shp = "data/05_profile_lines" + ".shp"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos los canales principales
    main_heads = p.heads_from_points(dem, main_ch)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Abrimos shapefile de poligonos con las cuencas
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)

    perfiles = []
    for feat in layer:
        # Get the basin and the heads inside
        basin_geom = feat.GetGeometryRef()
        basin_heads = p.heads_inside_basin(heads, basin_geom)

        # Get profiles
        pps = p.get_profiles(fac, dem, basin_heads, basin_geom, tributaries=True)
        perfiles.extend(pps)

    p.profiles_to_shp(out_point_shp, perfiles)
    p.profiles_to_shp(out_line_shp, perfiles, dl)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


test01()
