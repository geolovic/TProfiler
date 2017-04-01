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
print "Tests for profiler.get_profiles()"

    
def test01():
    """
    Test for get_profiles() function
    """
    inicio = time.time()
    print "=" * 40
    print "Test 01 para get_profiles() function"
    print "Test in progress..."
    
    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    umbral = 1000
    units = "CELL"
    out_shp = "data/03_profiles_test01.shp"

    cabeceras = p.get_heads(fac, dem, umbral, units=units)
    perfiles = p.get_profiles(fac, dem, cabeceras, tributaries=True)
    p.profiles_to_shp(out_shp, perfiles)
    
    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40


def test02():
    """
    Test for get_profiles() function
    Testeando una sola cuenca
    """
    inicio = time.time()
    print "=" * 40
    print "Test 02 para get_profiles() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    umbral = 1000
    units = "CELL"
    out_shp = "data/03_profiles_test02.shp"

    # Get input basin geometry
    input_basin_shp = "data/cuenca_darro.shp"
    dataset = ogr.Open(input_basin_shp)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(0)
    basin = feat.GetGeometryRef()

    # Get heads and profiles
    cabeceras = p.get_heads(fac, dem, umbral, units=units)
    cabeceras = p.heads_inside_basin(cabeceras, basin)
    perfiles = p.get_profiles(fac, dem, cabeceras, basin, tributaries=True)
    p.profiles_to_shp(out_shp, perfiles)

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40


def test03():
    """
    Test for get_profiles() function
    Testeando todas las cuencas de un shapefile
    """
    inicio = time.time()
    print "=" * 40
    print "Test 03 para get_profiles() function"
    print "Test in progress..."

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    umbral = 1000
    units = "CELL"
    out_shp = "data/03_profiles_test03.shp"

    # Get input basin geometry
    input_basin_shp = "data/cuencas.shp"
    dataset = ogr.Open(input_basin_shp)
    layer = dataset.GetLayer(0)

    # Get profiles
    perfiles = []
    for feat in layer:
        basin = feat.GetGeometryRef()
        cabeceras = p.get_heads(fac, dem, umbral, units=units)
        cabeceras = p.heads_inside_basin(cabeceras, basin)
        perfiles.extend(p.get_profiles(fac, dem, cabeceras, basin, tributaries=True))

    p.profiles_to_shp(out_shp, perfiles)

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40

test01()
test02()
test03()