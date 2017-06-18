# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
March, 2017
Testing suite for profiler.py
Last modified: 13 March 2017
"""

import time
import profiler as p
import ogr
import matplotlib.pyplot as plt
print "Tests for profiler.profiles_to_shp()"

    
def test01():
    """
    Test for profiles_to_shp() function (point shapefile)
    """
    inicio = time.time()
    print "=" * 40
    print "Test 01 para profiles_to_shp() function"
    print "Output point and line shapefiles"
    print "Test in progress..."
    
    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuencas.shp"
    main_ch = "data/main_heads.shp"
    umbral = 1000
    units = "CELL"
    out_point_shp = "data/04_profile_points" + ".shp"
    out_line_shp = "data/04_profile_lines" + ".shp"

    # Obtenemos todas las cabeceras del DEM y los canales principales
    heads = p.get_heads(fac, dem, umbral, units)
    main_heads = p.heads_from_points(dem, main_ch)

    # Abrimos shapefile de poligonos con las cuencas
    dataset = ogr.Open(basin)
    layer = dataset.GetLayer(0)

    perfiles = []
    for feat in layer:
        # Get the basin and the heads inside
        basin_geom = feat.GetGeometryRef()
        basin_heads = p.heads_inside_basin(heads, basin_geom)
        main_basin_heads = p.heads_inside_basin(basin_heads, basin_geom)
        main_basin_heads.reverse()
        for hh in main_basin_heads:
            basin_heads.insert(0, hh)

        # Get profiles
        pps = p.get_profiles(fac, dem, basin_heads, basin_geom, tributaries=True)
        perfiles.extend(pps)

    p.profiles_to_shp(out_point_shp, perfiles)
    p.profiles_to_shp(out_line_shp, perfiles, 250)

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40


test01()
