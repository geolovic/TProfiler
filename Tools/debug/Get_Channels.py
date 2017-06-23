# -*- coding: iso-8859-15 -*-
#
#  Get_Channels.py
#
#  Copyright (C) 2017  J. Vicente Perez, Universidad de Granada
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com

#  Version: 1.0
#  June 21, 2017

#  Last modified June 21, 2017

import ogr
import osr
import gdal
import numpy as np
from profiler import get_heads, heads_from_points, heads_inside_basin
from praster import open_raster, create_from_template


# DEBUG ARGUMENTS
# ===============

basedir = "../../test/data/"
dem = basedir + "darro25.tif"
fac = basedir + "darro25fac.tif"
out_channels = basedir + "out_main_ch.shp"
threshold = 1000
units = "CELL"
basin_shapefile = basedir + "cuencas.shp"
head_shapefile = basedir + "main_heads.shp"
id_field = "id"


# PROGRAM CODE
# =============

NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}


def main(dem, fac, out_channels, threshold, units, basin_shp, head_shapefile, id_field):

    # Create output shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(out_channels)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(open_raster(dem).proj)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Get all the heads
    if threshold:
        heads = get_heads(fac, dem, threshold, units)
    else:
        heads = np.array([], dtype="float32").reshape(0, 6)

    if head_shapefile:
        main_heads = heads_from_points(dem, head_shapefile, id_field)
        heads = np.append(main_heads, heads, axis=0)

    if heads.shape[0] == 0:
        return

    # Get all the channels
    channels = []
    if basin_shp:
        basin_dataset = driver.Open(basin_shp)
        basin_layer = basin_dataset.GetLayer(0)
        for feat in basin_layer:
            basin = feat.GetGeometryRef()
            basin_heads = heads_inside_basin(heads, basin)
            channels.extend(get_channels(fac, dem, basin_heads, basin))
    else:
        channels = get_channels(fac, dem, heads)

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
    facraster = open_raster(fac)
    demraster = open_raster(dem)

    # Auxiliar PRaster object to record proccesed cells
    aux_raster = create_from_template(dem, gdal.GDT_Byte, nodata=0)

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
                    chandata.append((next_point[0], next_point[1], z))
                    aux_raster.set_cell_value(next_pos, 1)
                    break

            # Se anaden a las listas
            chandata.append((next_point[0], next_point[1], z))

            # Se comprueba si esta "marcado";
            # si lo esta, coge el valor de chi0 y termina el bucle
            # sino, se marca (con un 1) y se coge el siguiente punto
            if aux_raster.get_cell_value(next_pos) == 1:
                break
            else:
                aux_raster.set_cell_value(next_pos, 1)

            next_pos = facraster.get_flow(next_pos)

        if len(chandata) > 5:
            out_channels.append((int(head[5]), np.array(chandata)))
        first_river = False

    return out_channels


main(dem, fac, out_channels, threshold=0, units=1000, basin_shp="", head_shapefile="", id_field="")
