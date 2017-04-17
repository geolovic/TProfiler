# -*- coding: iso-8859-15 -*-
#
#  Chi_Map_import.py
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
#  April 16, 2017

#  Last modified April 16, 2017

import ogr
import osr
import gdal
import numpy as np

##DEM=raster
##Flow_accumulation=raster
##Output_channels=output vector
##Use_threshold=boolean True
##Threshold=number 1000
##Units=selection CELLS;MAP
##Use_basins=boolean True
##Basin_shapefile=vector
##Use_heads=boolean True
##Heads_shapefile=vector
##Names_field=field Heads_shapefile

# Debug
DEM = "../../test/data/darro25.tif"
Flow_accumulation = "../../test/data/darro25fac.tif"
Output_channels = "../../test/data/out_main_ch.shp"
Use_threshold = True
Threshold = 1000
Units = "CELL"
Use_basins = True
Basin_shapefile = "../../test/data/cuencas.shp"
Use_heads = True
Heads_shapefile = "../../test/data/main_channels.shp"
Names_field = "name"
# End debug


NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}


class QgsIface:

    def setText(self, mensaje):
        print mensaje

    def setPercent(self, percent):
        print str(percent * 100) + "%"


def main(dem, fac, use_umbral, umbral, units, use_basins, cuencas_shp, use_heads, main_channels, name, out_channels):
    progress = QgsIface()

    # Get all heads
    progress.setText("Extracting heads ...")
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
    sp.ImportFromWkt(open_raster(dem).proj)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("Name", ogr.OFTString))

    progress.setText("Processing channels ...")
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
    facraster = open_raster(fac)
    demraster = open_raster(dem)

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
    demraster = open_raster(dem)
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
            name = str(feat.GetField(names_field))
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
    facraster = open_raster(fac)
    demraster = open_raster(dem)

    # Auxiliar PRaster object to record values of Chi
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

        if len(chandata) > 5:
            out_channels.append((head[5], np.array(chandata)))
        first_river = False
    return out_channels


def open_raster(raster_path):
    """
    This function open a raster and returns a pRaster instance

    :param raster_path: [str] Path to the raster to open
    :return: pRaster class instance
    """
    raster = gdal.Open(raster_path)
    if not raster:
        return
    array = raster.GetRasterBand(1).ReadAsArray()
    geot = raster.GetGeoTransform()
    proj = raster.GetProjection()
    nodata = raster.GetRasterBand(1).GetNoDataValue()

    return PRaster(array, geot, proj, nodata)


def create_raster(xsize, ysize, dtype=gdal.GDT_Int16, proj="", geot=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0), nodata=0.0):
    """
    This function creates a pRaster object "In Memory", to save it use the method Save(path)

    :param xsize: *int* -- Number of columns of the raster
    :param ysize: *int* -- Number of rows of the raster
    :param dtype: *gdal.GDT type* -- Data type of the new raster (Default = gdal.GDT_Int16)
    :param proj:  *str* -- Projection of the new raster in wkt (Default = "")
    :param geot:  *tuple* -- Geotransform matrix for the new raster (Default = (0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    :param nodata: *float* -- Nodata value for the new raster (Default = 0.0)
    :return: pRaster instance
    """

    # Creates an empty array and fill up with nodata values
    arrdata = np.empty((ysize, xsize)).astype(GTYPES[dtype])
    arrdata.fill(nodata)

    return PRaster(arrdata, proj, geot, nodata)


def create_from_template(template, dtype=None, nodata=None):
    """
    This function creates a raster with the same parameters than the template. The created raster is "In Memory", to
    save it use the method Save(path)

    :param template: *str* -- Path to the raster template
    :param dtype:  *gdal.GDT type* -- Data type for the new raster (Default = None -> Takes dtype from template)
    :param nodata: *float* -- NoData value for the new raster (Default = None -> Takes nodata from template)
    :return: pRaster instance
    """
    temp_raster = gdal.Open(template)
    temp_banda = temp_raster.GetRasterBand(1)
    geot = temp_raster.GetGeoTransform()
    proj = temp_raster.GetProjection()
    xsize = temp_banda.XSize
    ysize = temp_banda.YSize

    if dtype is None:
        dtype = temp_banda.DataType

    if nodata is None:
        nodata = temp_banda.GetNoDataValue()

    arrdata = np.empty((ysize, xsize)).astype(GTYPES[dtype])
    arrdata.fill(nodata)

    return PRaster(arrdata, geot, proj, nodata)


class PRaster:

    def __init__(self, array, geot=(0.0, 1.0, 0.0, 0.0, 0.0, -1.0), proj="", nodata=None):
        """
        Class to manipulate Raster objects

        :param array: *numpy.ndarray* -- Numpy array with raster data
        :param geot:  *tuple* -- Geotramsformation Matrix (upX, xcell, 0, upY, 0, ycell) (Default = (0,1,0,0,0,-1))
        :param proj:  *str* -- Projection of the new raster in wkt (Default = "")
        :param nodata: *float* -- Nodata value for the new raster (Default = None)
        """
        self.geot = geot
        self.proj = proj
        self.cellsize = geot[1]
        self.array = array
        self.nodata = nodata
        self.YSize = array.shape[0]
        self.XSize = array.shape[1]
        self.XMin = geot[0]
        self.YMin = geot[3] - self.cellsize * self.YSize
        self.XMax = geot[0] + self.cellsize * self.XSize
        self.YMax = geot[3]

    def get_cell_value(self, cell):
        """
        Get the raster value at a cell location

        :param cell: *tuple* -- (row, col) Cell location
        :return: Value of raster in the cell location
        """
        if cell[0] < 0 or cell[1] < 0:
            return None
        elif cell[0] >= self.YSize or cell[1] >= self.XSize:
            return None
        else:
            return self.array[cell[0], cell[1]]

    def get_xy_value(self, point):
        """
        Get the raster value at a point location

        :param point: *tuple* -- (X, Y) Point location
        :return: Value of raster in the point location
        """
        cell = self.xy_2_cell(point)
        return self.get_cell_value(cell)

    def xy_2_cell(self, point):
        """
        Get the cell position (row, col) of a point

        :param point: *tuple* -- (X, Y) Point location
        :return: Cell position (row, col)
        """
        row = int((self.YMax - point[1]) / self.cellsize)
        col = int((point[0] - self.XMin) / self.cellsize)
        return row, col

    def cell_2_xy(self, cell):
        """
        Get the XY position (X, Y) of a raster cell

        :param cell: *tuple* -- (row, col) Cell location
        :return: XY position (X, Y)
        """
        x = self.XMin + self.cellsize * cell[1] + self.cellsize / 2.
        y = self.YMax - self.cellsize * cell[0] - self.cellsize / 2.
        return x, y

    def set_cell_value(self, cell, value):
        """
        Set the raster value at a cell location

        :param cell: *tuple* -- (row, col) Cell location
        :param value: *number* -- New value for the pRaster at cell location
        """
        self.array[cell[0], cell[1]] = value

    def get_window(self, cell, ncells):
        """
        Get a ncells x ncells window around an specific cell. Function includes edge treatment.

        :param cell: *tuple* -- (row, col) Cell location
        :param ncells: *int* -- Number of cells around cell (Full window size will be 2*ncells + 1)
        :return: numpy.ndarray with the data of the raster within the window
        """
        row0 = cell[0] - ncells
        col0 = cell[1] - ncells

        nrows = ncells * 2 + 1
        ncols = ncells * 2 + 1

        # Check boundary conditions
        if row0 < 0:
            nrows += row0
            row0 = 0
        if col0 < 0:
            ncols += col0
            col0 = 0

        if row0 + nrows >= self.array.shape[0]:
            nrows = self.array.shape[0] - row0
        if col0 + ncols >= self.array.shape[1]:
            ncols = self.array.shape[1] - col0

            # Get raster data in the window and return an array
        window = self.array[row0:row0 + nrows, col0:col0 + ncols]

        if window.size > 0:
            return window
        else:
            return None

    def get_flow(self, cell):
        """
        Get the next cell where the flow goes (in case pRaster is a Flow Accumulation raster)

        :param cell: *tuple* -- (row, col) Cell location
        :return: cell (row, col) where the flow goes. It returns None at the edges of the raster
        """

        vec_adyacentes = [(-1, 0), (0, -1), (0, 1), (1, 0)]
        vec_diagonales = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

        # Suponemos que el valor máximo es el mismo
        cell_value = self.get_cell_value(cell)

        max_value = cell_value
        max_pos = cell

        # La celda a la que va el flujo no tiene porque ser la de mayor valor de flow accumulation
        # En el caso de que el flujo haga una L la máxima es la diagonal, pero el flujo va a la adyacente
        # Por ello primero se comprueban las celdas adyacentes y luego las diagonales

        for n in vec_adyacentes:
            row = cell[0] + n[0]
            col = cell[1] + n[1]
            if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                # A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                if self.get_cell_value((row, col)) > max_value:
                    max_value = self.get_cell_value((row, col))
                    max_pos = (row, col)

        if cell_value == max_value:
            # Si no hay ninguna celda adyacente con un valor mayor de f
            for n in vec_diagonales:
                row = cell[0] + n[0]
                col = cell[1] + n[1]
                if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                    # A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                    if self.get_cell_value((row, col)) > max_value and max_value != self.nodata:
                        max_value = self.get_cell_value((row, col))
                        max_pos = (row, col)

        if max_value == cell_value or max_value == self.nodata:
            return None
        else:
            return max_pos

    def save_raster(self, path):
        """
        Saves the pRaster in the disk

        :param path: *str* -- Path where new raster will be saved
        """

        if self.array.dtype not in NTYPES.keys():
            return
        else:
            tipo = NTYPES[str(self.array.dtype)]

        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self.XSize, self.YSize, 1, tipo)
        raster.SetGeoTransform(self.geot)
        raster.set_projection(self.proj)
        if self.nodata:
            raster.GetRasterBand(1).SetNoDataValue(self.nodata)
        raster.GetRasterBand(1).WriteArray(self.array)


main(DEM, Flow_accumulation, Use_threshold, Threshold, Units, Use_basins, Basin_shapefile,
     Use_heads, Heads_shapefile, Names_field, Output_channels)
