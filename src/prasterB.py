# -*- coding: utf-8 -*-
#
#  prasterB.py
#
#  Copyright (C) 2016  J. Vicente Perez, Universidad de Granada
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

#  Version: 3.0
#  November 23, 2016

#  Last modified December 6, 2016

import tempfile
import gdal
import struct

DATA_TYPES = {1:"h", 2:"H", 3:"h", 4:"I", 5:"i", 6:"f",
              7:"D", 8:"h", 9:"i", 10:"f", 11:"D"}

def Open(raster_path):
    """
    This function open a raster and returns a pRaster instance

    :param raster_path: [str] Path to the raster to open
    :return: pRaster class instance
    """
    return pRaster(gdal.Open(raster_path, gdal.GA_Update))

def Create(xsize, ysize, dtype = gdal.GDT_Int16, proj = "", geot = (0.0, 1.0, 0.0, 0.0, 0.0, 1.0), nodata = 0.0):
    """
    This function creates a pRaster object in a temporary folder, to save it in other place, use the method Save(path)

    :param xsize: *int* -- Number of columns of the raster
    :param ysize: *int* -- Number of rows of the raster
    :param dtype: *gdal.GDT type* -- Data type of the new raster (Default = gdal.GDT_Int16)
    :param proj:  *str* -- Projection of the new raster in wkt (Default = "")
    :param geot:  *tuple* -- Geotransform matrix for the new raster (Default = (0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    :param nodata: *float* -- Nodata value for the new raster (Default = 0.0)
    :return: pRaster instance
    """

    tmp_dir = tempfile.mkdtemp(prefix='prasTMP_')
    tmp_name = next(tempfile._get_candidate_names())
    path = tmp_dir + "/" + tmp_name + ".tif"

    driver = gdal.GetDriverByName("GTiff")
    raster = driver.Create(path, xsize, ysize, 1, dtype)
    raster.SetGeoTransform(geot)
    raster.set_projection(proj)
    raster.GetRasterBand(1).SetNoDataValue(nodata)
    
    return pRaster(raster)

def CreateFromTemplate(template, dtype = None, nodata = None):
    """
    This function creates a raster with the same parameters than the template in a temporary directory, to
    save it use the method Save(path)

    :param template: *str* -- Path to the raster template
    :param dtype:  *gdal.GDT type* -- Data type for the new raster (Default = None -> Takes dtype from template)
    :param nodata: *float* -- NoData value for the new raster (Default = None -> Takes nodata from template)
    :return: pRaster instance
    """

    driver = gdal.GetDriverByName("GTiff")
    temp_raster = gdal.Open(template)
    temp_banda = temp_raster.GetRasterBand(1)

    if dtype is None:
        dtype = temp_banda.DataType

    if nodata is None:
        nodata = temp_banda.GetNoDataValue()

    tmp_dir = tempfile.mkdtemp(prefix='prasTMP_')
    tmp_name = next(tempfile._get_candidate_names())
    path = tmp_dir + "/" + tmp_name + ".tif"

    xsize = temp_banda.XSize
    ysize = temp_banda.YSize
    raster = driver.Create(path, xsize, ysize, 1, dtype)
    raster.SetGeoTransform(temp_raster.GetGeoTransform())
    raster.set_projection(temp_raster.get_projection())
    raster.GetRasterBand(1).SetNoDataValue(nodata)

    return pRaster(raster)

class pRaster:
    """
    Class to manipulate Raster objects

    :param in_raster: *gdal.raster* -- gdal.raster object
    """    
    
    def __init__(self, in_raster):
        self.raster = in_raster
        geot = self.raster.GetGeoTransform()
        self.banda = self.raster.GetRasterBand(1)
        self.cellsize = geot[1]
        self.nodata = self.banda.GetNoDataValue()
        self.YSize = self.banda.YSize
        self.XSize = self.banda.XSize       
        self.XMin = geot[0]
        self.YMin = geot[3] - self.cellsize * self.YSize
        self.XMax = geot[0] + self.cellsize * self.XSize
        self.YMax = geot[3]
        
    def GetCellValue(self, cell):
        """
        Get the raster value at a cell location

        :param cell: *tuple* -- (row, col) Cell location
        :return: Value of raster in the cell location
        """

        if cell[0] < 0 or cell[1] < 0:
            return None
        elif cell[0] >= self.YSize or cell[1]>= self.XSize:
            return None
        else:
            data = self.banda.ReadRaster(cell[1], cell[0], 1, 1, 1, 1, self.banda.DataType)
            return struct.unpack("<" + DATA_TYPES[self.banda.DataType], data)[0]

    def GetXYValue(self, point):
        """
        Get the raster value at a point location

        :param point: *tuple* -- (X, Y) Point location
        :return: Value of raster in the point location
        """

        cell = self.XY2Cell(point)
        return self.GetCellValue(cell)
 
    def XY2Cell(self, point):
        """
        Get the cell position (row, col) of a point

        :param point: *tuple* -- (X, Y) Point location
        :return: Cell position (row, col)
        """

        row = int((self.YMax - point[1]) / self.cellsize)
        col = int((point[0] - self.XMin) / self.cellsize)
        return (row, col)
    
    def Cell2XY(self, cell):
        """
        Get the XY position (X, Y) of a raster cell

        :param cell: *tuple* -- (row, col) Cell location
        :return: XY position (X, Y)
        """

        x = self.XMin + (self.cellsize * cell[1]) + (self.cellsize / 2)
        y = self.YMax - (self.cellsize * cell[0]) - (self.cellsize / 2)
        return (x, y)

    def SetCellValue(self, cell, value):
        """
        Set the raster value at a cell location

        :param cell: *tuple* -- (row, col) Cell location
        :param value: *number* -- New value for the pRaster at cell location
        """

        data = struct.pack("<" + DATA_TYPES[self.banda.DataType], value)
        self.banda.WriteRaster(cell[1], cell[0], 1, 1, data)

    def GetWindow(self, cell, ncells):
        """
        Get a ncells x ncells window around an specific cell. Function includes edge treatment.

        :param cell: *tuple* -- (row, col) Cell location
        :param ncells: *int* -- Number of cells around cell (Full window size will be 2*ncells + 1)
        :return: numpy.ndarray with the data of the raster within the window
        """

        xoff = cell[1] - ncells
        yoff = cell[0] - ncells        
        xcells = ncells * 2 + 1
        ycells = ncells * 2 + 1
        
        #Check boundary conditions
        if xoff < 0:
            xcells = xcells + xoff
            xoff = 0   
        if yoff < 0:
            ycells = ycells + yoff
            yoff = 0
            
        if xoff + xcells >= self.banda.XSize:
            xcells = self.banda.XSize - xoff
        if yoff + ycells >= self.banda.YSize:
            ycells = self.banda.YSize - yoff    
        
        #Get raster data in the window and return an array
        data = self.banda.ReadAsArray(xoff, yoff, xcells, ycells)
        
        return data
        
    def GetFlow(self, cell):
        """
        Get the next cell where the flow goes (in case pRaster is a Flow Accumulation raster)

        :param cell: *tuple* -- (row, col) Cell location
        :return: cell (row, col) where the flow goes. It returns None at the edges of the raster
        """

        vec_adyacentes = [(-1,0),(0,-1),(0,1),(1,0)]
        vec_diagonales = [(-1,-1),(-1,1),(1,-1),(1,1)]
        
        #Suponemos que el valor maximo es el mismo
        cell_value = self.GetCellValue(cell)
        max_value = cell_value
        max_pos = cell
    
        #La celda a la que va el flujo no tiene porque ser la de mayor valor de flow accumulation
        #En el caso de que el flujo haga una L la máxima es la diagonal, pero el flujo va a la adyacente
        #Por ello primero se comprueban las celdas adyacentes y luego las diagonales
    
        for n in vec_adyacentes:
            row = cell[0] + n[0]
            col = cell[1] + n[1]
            if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                #A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                if self.GetCellValue((row, col)) > max_value:
                    max_value = self.GetCellValue((row, col))
                    max_pos = (row, col)
        
        if cell_value == max_value:
            #Si no hay ninguna celda adyacente con un valor mayor de f
            for n in vec_diagonales:
                row = cell[0] + n[0]
                col = cell[1] + n[1]
                if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                    #A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                    if self.GetCellValue((row, col)) > max_value and max_value != self.nodata:
                        max_value = self.GetCellValue((row, col))
                        max_pos = (row, col)            
        
        if max_value == cell_value or max_value == self.nodata:
            return None
        else:
            return max_pos

    def Save(self, path):
        """
        Saves the pRaster in the disk

        :param path: *str* -- Path where new raster will be saved
        """
        driver = gdal.GetDriverByName("GTiff")
        dst_ras = driver.CreateCopy(path, self.raster, 0)