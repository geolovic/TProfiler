# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
July, 2017

Testing suite for get_channels.py
Last modified: 13 July 2017
"""
import ogr
import osr
import gdal
import os
import profiler as p
import time
import numpy as np


def test_01():
    """
    Test for profiles_to_shp() function
    Extracting all the channels in a DEM
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para get_channels() function")
    print("Extrating all the channels in a DEM")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    umbral = 1000
    units = "CELL"
    out_channels = "data/out/07_out_channels_01.shp"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todos los canales
    channels = p.get_channels(fac, dem, heads)

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


def test_02():
    """
    Test for profiles_to_shp() function
    Extracting all the channels in a DEM using main channels
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para get_channels() function")
    print("Extrating all the channels in a DEM but using main heads")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    umbral = 1000
    units = "CELL"
    main_channels = "data/in/main_heads.shp"
    id_field = "id"
    out_channels = "data/out/07_out_channels_02.shp"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_channels, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos todos los canales
    channels = p.get_channels(fac, dem, heads)

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


def test_03():
    """
    Test for profiles_to_shp() function
    Extracting all the channels in a DEM inside basins
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 03 para get_channels() function")
    print("Extrating all the channels in a DEM inside basins")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    umbral = 1000
    units = "CELL"
    basin_shp = "data/in/cuencas.shp"
    out_channels = "data/out/07_out_channels_03.shp"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos los diferentes poligonos de las cuencas
    dataset = ogr.Open(basin_shp)
    layer = dataset.GetLayer(0)

    # Otenemos los canales dentro de cada cuenca
    channels = []
    for feat in layer:
        basin_geom = feat.GetGeometryRef()
        heads_inside = p.heads_inside_basin(heads, basin_geom)
        channels.extend(p.get_channels(fac, dem, heads_inside, basin_geom))

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)

    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


def test_04():
    """
    Test for profiles_to_shp() function
    Extracting all the channels in a DEM inside basins and using main channels
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 04 para get_channels() function")
    print("Extrating all the channels in a DEM inside basins and using main channels")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    umbral = 1000
    units = "CELL"
    basin_shp = "data/in/cuencas.shp"
    main_channels = "data/in/main_heads.shp"
    id_field = "id"
    out_channels = "data/out/07_out_channels_04.shp"

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_channels, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos los diferentes poligonos de las cuencas
    dataset = ogr.Open(basin_shp)
    layer = dataset.GetLayer(0)

    # Otenemos los canales dentro de cada cuenca
    channels = []
    for feat in layer:
        basin_geom = feat.GetGeometryRef()
        heads_inside = p.heads_inside_basin(heads, basin_geom)
        channels.extend(p.get_channels(fac, dem, heads_inside, basin_geom))

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)

    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


def test_05():
    """
    Test for profiles_to_shp() function
    Extracting only selected channels for the entire DEM
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 05 para get_channels() function")
    print("Extrating only selected channels for the entire DEM")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    main_channels = "data/in/main_heads.shp"
    id_field = "id"
    out_channels = "data/out/07_out_channels_05.shp"

    # Obtenemos las cabeceras de los rios principales
    heads = p.heads_from_points(dem, main_channels, id_field=id_field)

    # Otenemos los canales de los rios principales
    channels = p.get_channels(fac, dem, heads)

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)

    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


def test_06():
    """
    Test for profiles_to_shp() function
    Extracting only selected channels inside basins
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 06 para get_channels() function")
    print("Extracting only selected channels inside basins")
    print("Test in progress...")

    # Test parameters
    dem = "data/in/darro25.tif"
    fac = "data/in/darro25fac.tif"
    basin_shp = "data/in/cuencas.shp"
    main_channels = "data/in/main_heads.shp"
    id_field = "id"
    out_channels = "data/out/07_out_channels_06.shp"

    # Obtenemos las cabeceras de los rios principales
    heads = p.heads_from_points(dem, main_channels, id_field=id_field)

    # Obtenemos los diferentes poligonos de las cuencas
    dataset = ogr.Open(basin_shp)
    layer = dataset.GetLayer(0)

    # Otenemos los canales dentro de cada cuenca
    channels = []
    for feat in layer:
        basin_geom = feat.GetGeometryRef()
        heads_inside = p.heads_inside_basin(heads, basin_geom)
        channels.extend(p.get_channels(fac, dem, heads_inside, basin_geom))

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_channels):
        dataset = driver.Open(out_channels, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_channels)

    sp = osr.SpatialReference()
    sp.ImportFromWkt(proj_wkt)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Guardamos channels en shapefile
    for ch in channels:
        name = ch[0]
        chandata = ch[1]
        geom = ogr.Geometry(ogr.wkbLineString)
        for row in chandata:
            geom.AddPoint(float(row[0]), float(row[1]))
        # Create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Add geometry to feature and feature to layer
        feature.SetGeometry(geom)
        feature.SetField("id", name)
        layer.CreateFeature(feature)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Output channels guardados en {0}".format(out_channels))
    print("=" * 40)


test_01()
test_02()
test_03()
test_04()
test_05()
test_06()
