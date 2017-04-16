# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
March, 2017
Testing suite for profiler.py
Last modified: 13 March 2017
"""

import praster as p
import ogr
import osr
import gdal
import numpy as np
print "Tests for profiler.get_channels"

DEM = r"data/darro25.tif"
Flow_accumulation = "data/darro25fac.tif"
Output_channels = "data/out_main_ch.shp"
Use_threshold = True
Threshold = 1000
Units = "CELL"
Use_basins = True
Basin_shapefile = "data/cuencas.shp"
Use_heads = True
Heads_shapefile = "data/main_channels.shp"
Names_field = "name"


def get_heads(fac, dem, umbral, units="CELL"):
    """
    Extracts heads from flow accumulation raster based in a threshold

    :param fac: *str* -- Path to the flow acculumation raster
    :param dem: *str* -- Path to the DEM
    :param umbral: *float* -- Threshold for channel initiation (i.e. for head definition)
    :param units: *str -- Units of threshold ("MAP" / "CELL" ) (Default="MAP")
    :return: *list* -- List of tuples (row, col, X, Y, Z, id)
    """

    # Open fac and dem rasters
    facraster = p.open_raster(fac)
    demraster = p.open_raster(dem)

    if units == "MAP":
        umbral = int(umbral / facraster.cellsize ** 2)
    else:
        umbral = int(umbral)

    heads = []

    # A cell will be a true head if its area >= threshold and it hasn't got
    # any neighbouring cell with a lower area than the position but higher than threshold
    selcells = np.where(facraster.array >= umbral)
    if selcells[0].size > 0:
        arr = np.array(selcells).T
        for cell in arr:
            row = int(cell[0])
            col = int(cell[1])
            area = facraster.get_cell_value((row, col))
            elev = demraster.get_cell_value((row, col))
            point = facraster.cell_2_xy((row, col))
            if area == umbral:
                heads.append((row, col, point[0], point[1], elev))
            elif (area >= umbral) and (area < umbral * 3):
                vecinos = facraster.get_window((row, col), 1)
                if len(np.where((vecinos >= umbral) & (vecinos < area))[0]) == 0:
                    heads.append((row, col, point[0], point[1], elev))

    # If there is any heads, return an empty list
    if len(heads) == 0:
        output_heads = heads
    else:
        # Sort heads by elevation
        heads = np.array(heads)
        ind = heads[:, 4].argsort()
        ind = ind[::-1]
        heads = heads[ind]

        # Add a column id (string) and return heads
        head_id = 0
        output_heads = []
        for hh in heads:
            output_heads.append((int(hh[0]), int(hh[1]), float(hh[2]), float(hh[3]), float(hh[4]), str(head_id)))
            head_id += 1

    return output_heads


def heads_from_points(dem, point_shp, names_field=""):
    """
    This function creates a list of tuples (row, col, X, Y, Z, "id") from a point shapefile

    :param dem: *str* -- Path to the DEM
    :param point_shp: *str* -- Path to the point shapefile
    :param names_field: *str* -- String with the field with profile names
    :return: *list* -- List of tuples (row, col, X, Y, Z, 'id')
    """

    heads = []
    demraster = p.open_raster(dem)
    dataset = ogr.Open(point_shp)
    layer = dataset.GetLayer(0)
    n = 0
    for feat in layer:
        geom = feat.GetGeometryRef()
        punto = (geom.GetX(), geom.GetY())
        cell = demraster.xy_2_cell(punto)
        elev = demraster.get_cell_value(cell)
        layerdef = layer.GetLayerDefn()
        fields = [layerdef.GetFieldDefn(idx).GetName() for idx in range(layerdef.GetFieldCount())]
        if names_field in fields:
            name = str(feat[names_field])
        else:
            name = str(n)
        n += 1
        heads.append([cell[0], cell[1], punto[0], punto[1], elev, name])

    return heads


def heads_inside_basin(heads, basin):
    """
    This function extracts heads that pass a threshold from a flow accumulation raster inside a determined basin

    :param heads: *list* -- List of tuples (row, col, X, Y, elev, id)
    :param basin: *ogr.Geometry(ogr.wkbLinearRing)* -- Polygon defining the basin
    :return: *list* -- List of tuples (row, col, X, Y, "id")
    """

    if len(heads) == 0:
        return heads

    basin_heads = []
    # Take only heads that are within basin polygon
    for head in heads:
        pto = ogr.Geometry(ogr.wkbPoint)
        pto.AddPoint(head[2], head[3])
        if basin.Contains(pto):
            basin_heads.append(head)

    return basin_heads


def get_channels(fac, dem, heads, basin=None):
    """
    This fucntion extracts profiles for all heads in heads list. It will complete river profiles until the edge of
    the dem or until the basin limits (if basin is specified).
    :param fac: PRaster Raster with the flow accumulation
    :param dem: PRaster Raster with the Digital Elevation Model
    :param heads: List of tuples (row, col, X, Y, Z, "id")
    :param basin: ogr.Geometry Polygon that represents one basin
    :return: list of arrays representing channels (3 columns: x, y, z)
    """

    # Get facraster and demrasters as pRaster objects
    facraster = p.open_raster(fac)
    demraster = p.open_raster(dem)

    # Auxiliar PRaster object to record values of Chi
    aux_raster = p.create_from_template(dem, gdal.GDT_Byte, nodata=0)

    # If all the heads come from the same basin, all of them will flow into first river
    # there is no need to check if all vertexes are inside the basin, from the second head all will be inside
    first_river = True

    out_channels = []
    # Process all the heads
    for head in heads:
        # List to store xyz tuples with channel data
        chandata = []

        # Calculate data for the first point
        pos = (int(head[0]), int(head[1]))
        point = demraster.cell_2_xy(pos)
        z = demraster.get_cell_value(pos)

        # Add first point to profile data (xyz)
        chandata.append((point[0], point[1], z))

        # 'Mark' aux raster to indicate that the pixel was processed
        aux_raster.set_cell_value(pos, 1)

        # Obtain the next point following the flow direction
        next_pos = facraster.get_flow(pos)

        # Start bucle to obtain vertexes following river's flow
        while next_pos:

            # Data of the new point
            next_point = demraster.cell_2_xy(next_pos)
            z = demraster.get_cell_value(next_pos)

            # Check if point is still inside basin (if specified and it is not the first head)
            if basin and first_river:
                pto = ogr.Geometry(ogr.wkbPoint)
                pto.AddPoint(next_point[0], next_point[1])
                if not basin.Contains(pto):
                    break

            # Se anaden a las listas
            chandata.append((next_point[0], next_point[1], z))

            # Se comprueba si esta "marcado";
            # si lo esta, coge el valor de chi0 y termina el bucle
            # sino, se marca (con un 1) y se coge el siguiente punto
            if aux_raster.get_cell_value(next_pos) == 1:
                break
            else:
                aux_raster.set_cell_value(pos, 1)
                pos = tuple(next_pos)
            next_pos = facraster.get_flow(next_pos)

        out_channels.append((head[5], np.array(chandata)))
        first_river = False
    return out_channels


def main(dem, fac, use_umbral, umbral, units, use_basins, cuencas_shp, use_heads, main_channels, name, out_channels):

    # Get all heads
    if use_umbral:
        heads = get_heads(fac, dem, umbral, units)
    else:
        heads = []

    if use_heads:
        main_heads = heads_from_points(dem, main_channels, name)
        for head in main_heads:
            heads.insert(0, head)

    # Get channels (each channel is a numpy array with 3 columns: x, y, z)
    # Open basin shapefile if specified
    if use_basins:
        channels = []
        dataset = ogr.Open(cuencas_shp)
        layer = dataset.GetLayer(0)
        for feat in layer:
            basin_geom = feat.GetGeometryRef()
            basin_heads = heads_inside_basin(heads, basin_geom)
            channels.extend(get_channels(fac, dem, basin_heads, basin_geom))
    else:
        channels = get_channels(fac, dem, heads)

    # Creamos output
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(out_channels)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(p.open_raster(dem).proj)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("Name", ogr.OFTString))

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
        feature.SetField("Name", name)
        layer.CreateFeature(feature)

main(DEM, Flow_accumulation, Use_threshold, Threshold, Units, Use_basins, Basin_shapefile,
     Use_heads, Heads_shapefile, Names_field, Output_channels)
