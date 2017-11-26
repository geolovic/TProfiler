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
#  July 14, 2017

#  Last modified 26 November, 2017

import ogr
import osr
import gdal
import numpy as np
import os
import profiler as p


# PROGRAM CODE
# =============
def main(dem, fac, out_shp, threshold=0, units="CELL", head_shp="", id_field="", basin_shp=""):

    # Obtenemos todas las cabeceras del DEM
    if threshold:
        heads = p.get_heads(fac, dem, threshold, units)
    else:
        heads = np.array([], dtype="float32").reshape((0, 6))

    # Obtenemos todas las cabeceras de la capa de puntos
    if head_shp:
        main_heads = p.heads_from_points(dem, head_shp, id_field=id_field)
    else:
        main_heads = np.array([], dtype="float32").reshape((0, 6))

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    if heads.shape[0] == 0:
        return
    
    # Numeramos nuevamente las cabeceras
    ord_id = np.arange(len(heads)).astype("float32")
    heads[:, 5] = ord_id

    # Obtenemos los diferentes poligonos de las cuencas
    if basin_shp:
        dataset = ogr.Open(basin_shp)
        layer = dataset.GetLayer(0)

        # Otenemos los canales dentro de cada cuenca
        channels = []
        for feat in layer:
            basin_geom = feat.GetGeometryRef()
            heads_inside = p.heads_inside_basin(heads, basin_geom)
            # Numeramos nuevamente las cabeceras
            ord_id = np.arange(len(heads_inside)).astype("float32")
            heads_inside[:, 5] = ord_id
            channels.extend(p.get_channels(fac, dem, heads_inside, basin_geom))
    else:
        channels = p.get_channels(fac, dem, heads)

    # Creamos output shapefile
    proj_wkt = gdal.Open(dem).GetProjection()
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out_shp):
        dataset = driver.Open(out_shp, 1)
        for n in range(dataset.GetLayerCount()):
            dataset.DeleteLayer(n)
    else:
        dataset = driver.CreateDataSource(out_shp)

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