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
##Id_field=field Heads_shapefile


NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}


def main(dem, fac, use_threshold, threshold, units, use_basins, basins_shp, use_heads, heads_shapefile, id_field, out_channels):

    # Create output shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(out_channels)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(open_raster(dem).proj)
    layer = dataset.CreateLayer("channels", sp, ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    # Get all the heads
    if use_threshold:
        heads = get_heads(fac, dem, threshold, units)
    else:
        heads = np.array([], dtype="float32").reshape(0, 6)

    if use_heads:
        main_heads = heads_from_points(dem, heads_shapefile, id_field)
        heads = np.append(main_heads, heads, axis=0)

    if heads.shape[0] == 0:
        return

    # Get all the channels
    channels = []
    if use_basins:
        basin_dataset = driver.Open(basins_shp)
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
            out_channels.append((int(head[5]), np.array(chandata)))
        first_river = False
    return out_channels


# IMPORTED FUNCTIONS (18 June 2017)
def get_heads(fac, dem, umbral, units="CELL"):
    """
    Extracts heads from flow accumulation raster based in a threshold

    :param fac: *str* -- Path to the flow acculumation raster
    :param dem: *str* -- Path to the DEM
    :param umbral: *float* -- Threshold for channel initiation (i.e. for head definition)
    :param units: *str -- Units of threshold ("MAP" / "CELL" ) (Default="MAP")
    :return: *numpy.array* -- Numpy.array with 6 columns [row, col, X, Y, Z, hid]
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

    # Si no hay ningun pixel que sea considerado una cabecera, devolvemos la lista vacia
    if len(heads) == 0:
        output_heads = heads
    else:
        nheads = len(heads)
        # Ordenamos las posiciones por su elevacion
        heads = np.array(heads, dtype="float32")
        ind = heads[:, 4].argsort()
        ind = ind[::-1]
        heads = heads[ind]

        # Anadimos columna con ids y devolvemos lista con cabeceras
        hids = np.arange(nheads).astype("float32").reshape(nheads, 1)
        output_heads = np.append(heads, hids, axis=1)

    return output_heads


def heads_from_points(dem, point_shp, id_field=""):
    """
    This function creates a list of tuples (row, col, X, Y, Z, "id") from a point shapefile

    :param dem: *str* -- Path to the DEM
    :param point_shp: *str* -- Path to the point shapefile
    :param id_field: *str* -- String with the id field (heads will be ordered by the values in this field)
    :return: *numpy.array* -- Numpy.array with 6 columns [row, col, X, Y, Z, hid]
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
        if id_field in fields:
            hid = str(feat[id_field])
        else:
            hid = n + 1
        n += 1
        heads.append([cell[0], cell[1], punto[0], punto[1], elev, hid])

    output_heads = np.array(heads, dtype="float32")
    ind = output_heads[:, 5].argsort()
    output_heads = output_heads[ind]

    return output_heads


def heads_inside_basin(heads, basin, first=0):
    """
    This function extracts heads that pass a threshold from a flow accumulation raster inside a determined basin

    :param heads: *numpy.array* -- Numpy.array with 6 columns [row, col, X, Y, Z, hid]
    :param basin: *ogr.Geometry(ogr.wkbLinearRing)* -- Polygon defining the basin
    :param first: *int* -- Identifier for the first head, all the other will have consecutives ids
    :return: *numpy.array* -- Numpy.array with 6 columns [row, col, X, Y, Z, hid]
    """

    if heads.shape[0] == 0:
        return heads

    basin_heads = []
    # Take only heads that are within basin polygon
    for head in heads:
        pto = ogr.Geometry(ogr.wkbPoint)
        pto.AddPoint(float(head[2]), float(head[3]))
        if basin.Contains(pto):
            basin_heads.append(head)

    output_heads = np.array(basin_heads, dtype="float32")
    output_heads[:, 5] = np.arange(output_heads.shape[0]).astype("float32") + float(first)
    return output_heads


class TProfile:
    """
    Properties:
    ============================
    self.dem_res :: *float*
      Dem resolution of the Digital elevation model used to retreive area-elevation data
    self.rid :: *int*
      Indentifier of the profile
    self.thetaref :: *float*
      Value of m/n used in area-slope calculations
    self.reg_points :: *int*
      Number of points used initially in slope and ksn regressions
    self.n_points :: *int*
      Number of vertexes of the profile
    self._data :: *numpy.array*
      11-column numpy.array with profile data

    =======   ==============================================
    Column    Description
    =======   ==============================================
    c0        X Coordinates of river profile vertex
    c1        Y Coordinates of river profile vertex
    c2        Z Elevation of the river profile vertex
    c3        L Distance to river head
    c4        A Drainage area to the vertex
    c5        Chi (Integral mode)
    c6        Slope of river profile in each vertex
    c7        ksn values for each vertex (calculated by linear regression)
    c8        Quality slope, correlation coefficient (r^2) of slope regressions
    c9        Quality ksn, correlation coefficient (r^2) of ksn regressions
    c10       Raw Z Elevation of the river profile vertex (used to reset the profile)
    =======   ==============================================
    """

    def __init__(self, pf_data, dem_res=0, rid=0, thetaref=0.45, chi0=0, reg_points=4, srs="", mouthdist=0, smooth=0):
        """
        Class that defines a river profile with morphometry capabilities.

        :param pf_data: *numpy array* - Array with input values
        :param dem_res: *float* - Resolution of the DEM used to extract profile features
        :param rid: *int* - Profile Identifier
        :param thetaref: *float* - Thetaref (m/n) value used to calculate Chi and Ksn indexes
        :param chi0: *float* - Value of chi index for first point (for tributaries)
        :param reg_points: *int* - Number of points (at each side) to calculate initial slope and ksn for each vertex
        :param srs: *str* - Spatial Reference system expresed as well knwon text (wkt)
        :param mouthdist: *float* - Distance from profile to the river mouth (for tributaries)
        :param smooth: *float* - Initial distance to smooth elevations (before to calculate slopes and ksn)

        pf_data param: (numpy.array) with at least 5 columns:

        =======   ==============================================
        Column    Description
        =======   ==============================================
        c0        X Coordinates of river profile vertex
        c1        Y Coordinates of river profile vertex
        c2        Z Elevation of the river profile vertex
        c3        L Distance to head (or to the first vertex)
        c4        A Drainage area of the vertex (in square meters!)
        =======   ==============================================
        """

        # Set profile properties
        self._srs = srs  # EPSG Code of the Spatial Reference
        self._mouthdist = mouthdist
        self.dem_res = float(dem_res)
        self.rid = rid
        self.thetaref = abs(thetaref)
        self.reg_points = reg_points
        self.n_points = pf_data.shape[0]

        # Get profile data from pf_data array
        aux_values = np.empty((pf_data.shape[0], 6))
        aux_values.fill(np.nan)
        self._data = np.append(pf_data, aux_values, axis=1)

        # Set raw elevations
        self._data[:, 10] = np.copy(self._data[:, 2])

        # Smooth profile elevations before to calculate ksn and chi
        self.smooth(smooth)

        # Create slopes, chi and ksn values
        self.calculate_slope(reg_points)
        self.calculate_chi(chi0=chi0)
        self.calculate_ksn(reg_points)

    def get_projection(self):
        """
        Returns a string with the projection as wkt
        """
        return self._srs

    def set_projection(self, projection):
        """
        Set the projection of the profile

        :param projection: str - String with the projection in wkt
        :return: None
        """
        self._srs = projection

    def length(self):
        """
        Returns the total length of the profile
        """
        return self._data[-1, 3]

    def get_x(self, head=True):
        """
        Returns x coordinates for all vertices

        :param head: boolean - Specifies if x coordinates are returned from head (True) or mouth (False)
        :return: numpy.array wiht x values for all vertices
        """
        if head:
            return np.copy(self._data[:, 0])
        else:
            return np.copy(self._data[::-1, 0])

    def get_y(self, head=True):
        """
        Returns y coordinates for all vertices

        :param head: boolean - Specifies if y coordinates are returned from head (True) or mouth (False)
        :return: numpy.array wiht y values for all vertices
        """
        if head:
            return np.copy(self._data[:, 1])
        else:
            return np.copy(self._data[::-1, 1])

    def get_z(self, head=True, relative=False):
        """
        Returns elevation values for all vertices

        :param head: boolean - Specifies if elevations are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if elevations are relative (min elevation = 0) or not
        :return: numpy.array wiht elevation values for all vertices
        """
        z_values = np.copy(self._data[:, 2])
        if relative:
            z_min = z_values[-1]
            z_values -= z_min

        if head:
            return z_values
        else:
            return z_values[::-1]

    def get_raw_z(self, head=True, relative=False):
        """
        Returns raw elevation values for all vertices

        :param head: boolean - Specifies if raw elevations are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if elevations are relative (min elevation = 0) or not
        :return: numpy.array wiht elevation values for all vertices
        """
        raw_z = np.copy(self._data[:, 10])
        if relative:
            z_min = raw_z[-1]
            raw_z -= z_min

        if head:
            return raw_z
        else:
            return raw_z[::-1]

    def get_l(self, head=True):
        """
        Returns a numpy.array with distances for all profile vertices

        :param head: boolean - Specifies if distances are returned from head (True) or mouth (False)
        If measured from mouth, a initial distance (self._mouthdist) will be added (to account for tributaries)
        :return: numpy.array with distances for all vertices (measured from head or mouth)
        """
        river_length = float(self._data[-1, 3])

        if head:
            li = np.copy(self._data[:, 3])
        else:
            li = river_length - self._data[:, 3] + self._mouthdist
            li = li[::-1]

        return li

    def get_area(self, head=True, cells=True):
        """
        Returns a numpy.array with drainage area values for all vertices

        :param head: boolean - Specifies if areas are returned from head (True) or mouth (False)
        :param cells: boolean - Specifies if areas are measured in cells (True) or in profile units (False)
        :return: numpy.array wiht area values for all vertices
        """
        areas = np.copy(self._data[:, 4])
        if not cells:
            areas *= self.dem_res ** 2

        if head:
            return areas
        else:
            return areas[::-1]

    def get_slope(self, threshold=0, head=True, lq=False):
        """
        Returns slopes calculated by linear regression

        :param threshold: *float* R^2 threshold. (Slopes with R^2 < threshold will be in lq_slopes array)
        :param head: boolean - Specifies if slopes are returned from head (True) or mouth (False)
        :param lq: boolean - Specifies lq_slopes will be returned or not
        :return: array with slopes (lq_slopes=False) or tuple of arrays (slopes, lq_slopes) (lq_slopes=True)
         slopes --> numpy.array of slopes with R^2 >= threshold (lq_slopes will receive a np.nan value)
         lq_slopes --> numpy.array of slopes with R^2 < threshold (slopes will receive a np.nan value)
        """
        slopes = []
        lq_slopes = []
        r2_values = self.get_slope_r2()
        for n in range(len(self._data)):
            if r2_values[n] >= threshold:
                slopes.append(self._data[n, 6])
                lq_slopes.append(np.nan)
            else:
                slopes.append(np.nan)
                lq_slopes.append(self._data[n, 6])

        slopes = np.array(slopes)
        lq_slopes = np.array(lq_slopes)

        if not head:
            slopes = slopes[::-1]
            lq_slopes = lq_slopes[::-1]

        if lq:
            return slopes, lq_slopes
        else:
            return slopes

    def get_ksn(self, threshold=0, head=True, lq=False):
        """
        Returns ksn values calculated by linear regression

        :param threshold: *float* R^2 threshold. (ksn with R^2 < threshold will be in lq_ksn array)
        :param head: boolean - Specifies if ksn are returned from head (True) or mouth (False)
        :param lq: boolean - Specifies lq_ksn will be returned or not
        :return: array with ksn (lq=False) or tuple of arrays (ksn, lq_ksn) (lq=True)
         ksn --> numpy.array of ksn values with R^2 >= threshold (lq_ksn will receive a np.nan value)
         lq_ksn --> numpy.array of ksn values with R^2 < threshold (ksn will receive a np.nan value)
        """
        ksn = []
        lq_ksn = []
        ksn_r2 = self.get_ksn_r2()
        for n in range(len(self._data)):
            if ksn_r2[n] >= threshold:
                ksn.append(self._data[n, 7])
                lq_ksn.append(np.nan)
            else:
                ksn.append(np.nan)
                lq_ksn.append(self._data[n, 7])

        ksn = np.array(ksn)
        lq_ksn = np.array(lq_ksn)

        if not head:
            ksn = ksn[::-1]
            lq_ksn = lq_ksn[::-1]

        if lq:
            return ksn, lq_ksn
        else:
            return ksn

    def get_slope_r2(self, head=True):
        """
        Returns slope R2 values from linear regressions for all vertices

        :param head: boolean - Specifies if R2 values are returned from head (True) or mouth (False)
        :return: numpy.array wiht slope R2 values for all vertices
        """
        if head:
            return np.copy(self._data[:, 8])
        else:
            return np.copy(self._data[::-1, 8])

    def get_ksn_r2(self, head=True):
        """
        Returns ksn R2 values from linear regressions for all vertices

        :param head: boolean - Specifies if R2 values are returned from head (True) or mouth (False)
        :return: numpy.array wiht ksn R2 values for all vertices
        """
        if head:
            return np.copy(self._data[:, 9])
        else:
            return np.copy(self._data[::-1, 9])

    def get_chi(self, head=True, relative=False):
        """
        Returns chi values for all vertices in ascending order.

        :param head: boolean - Specifies if chi values are returned from head (True) or mouth (False)
        :param relative: boolean - Specifies if chi values are relative (min chi = 0) or not
        :return: numpy.array wiht chi values for all vertices
        """
        chi_values = np.copy(self._data[:, 5])
        if relative:
            chi0 = chi_values[-1]
            chi_values -= chi0

        if head:
            return chi_values
        else:
            return chi_values[::-1]

    def smooth(self, window=0):
        """
        Smooths the elevations of the profile with a movil mean of window size. It also removes peaks and flat segments
        from the profile (to avoid problems when calculating slopes)

        :param window: Window size (in profile units) to smooth the elevations of the river profile
        :return: None
        """
        # Remove peaks and flat segments
        for n in range(len(self._data) - 1):
            if self._data[n + 1, 2] >= self._data[n, 2]:
                self._data[n + 1, 2] = float(self._data[n, 2]) - 0.001

        # Smooth elevations if window distance > 0
        if window > 0:
            n_cells = int(int((window / self.dem_res) + 0.5) / 2)
            for ind in range(len(self._data)):
                low = ind - n_cells
                high = ind + n_cells + 1
                if low < 0:
                    low = 0

                elevations = self._data[low:high, 10]
                self._data[ind, 2] = np.mean(elevations)

    def reset_elevations(self):
        """
        Reset smooth elevations. When reset, smooth elevations will equal to raw elevations
        """
        for n in range(len(self._data)):
            self._data[n, 2] = np.copy(self._data[n, 10])

    def calculate_chi(self, a0=1, chi0=0.0):
        """
        This function creates the chi data array. Chi data will be calculated for each vertex of the river profile and
        stored in the 7th column of self._data numpy array

        :param a0: *int* - Reference area to remove dimensionality of Chi index
        :param chi0: *float* - Initial Chi value in profile mouth. Needed to calculate chi values for tributaries
        """
        # Invert area array
        ai = self._data[::-1, 4] ** self.thetaref
        a0 = a0 ** self.thetaref
        chi = [chi0]
        for n in range(len(ai)):
            if n > 0:
                dx = self._data[n, 3] - self._data[n - 1, 3]
                chi.append(chi[n - 1] + (a0 * dx / ai[n]))

        self._data[:, 5] = chi[::-1]

    def calculate_slope(self, reg_points=4, raw_z=False):
        """
        This function calculates slopes for all vertexes by linear regression of distance-elevation data.
        Slopes are stored in column c5 of self._data. Together with slopes, R^2 are calculated (column  c6)

        :param reg_points: Number of profile points before and after each vertex to calculate slope
        :param raw_z: bool Specifies if raw_z values are taken from calculate slopes (True) or not (False)
        :return: None
        """
        li = self.get_l()
        if raw_z:
            zi = self.get_raw_z()
        else:
            zi = self.get_z()

        for n in range(self.n_points):
            low = n - reg_points
            high = n + reg_points

            if low < 0:
                low = 0

            sample_l = li[low:high + 1]
            sample_z = zi[low:high + 1]

            a = np.array([sample_l, np.ones(len(sample_l))]).T
            y = sample_z
            model, resid = np.linalg.lstsq(a, y)[:2]

            r2 = 1 - resid / (y.size * y.var())
            gradient = model[0]

            self._data[n, 8] = abs(r2)

            if abs(gradient) < 0.001:
                self._data[n, 6] = 0.001
            else:
                self._data[n, 6] = abs(gradient)

    def calculate_ksn(self, reg_points=4, raw_z=False):
        """
        This function calculates ksn for all vertexes by linear regression of chi-elevation data.

        :param reg_points: *int* - Number of profile points before and after each vertex to calculate ksn
        :param raw_z: bool Specifies if raw_z values are taken from calculate slopes (True) or not (False)
        :return: numpy.array with ksn values for all vertexes. If full is true, it returns a tuple of arrays (ksn, r^2)
        """
        ksn_values = []
        ksn_r2_values = []
        chi = self.get_chi(False)
        if raw_z:
            zi = self.get_raw_z(False)
        else:
            zi = self.get_z(False)

        for n in range(self.n_points):
            low = n - reg_points
            high = n + reg_points

            if low < 0:
                low = 0

            sample_chi = chi[low:high + 1]
            sample_z = zi[low:high + 1]

            poli, sce = np.polyfit(sample_chi, sample_z, deg=1, full=True)[:2]
            gradient = poli[0]
            r2 = 1 - sce / (sample_z.size * sample_z.var())

            ksn_r2_values.append(float(abs(r2)))

            if abs(gradient) < 0.0001:
                ksn_values.append(0.0001)
            else:
                ksn_values.append(abs(gradient))

        self._data[:, 7] = np.array(ksn_values)[::-1]
        self._data[:, 9] = np.array(ksn_r2_values)[::-1]

    def get_best_theta(self, a0=1, step=0.05):
        """
        Description
        ===========
        This function obtain the best m/n value for the profile following the approach
        proposed in Perron and Royden, 2013. This best m/n value will be the one that
        increases the linearity of the Chi-Elevation profile

        Parameters:
        ==============
        a0 :: *int (Default = 1)*
          Reference area value. By default set as 1 square meter

        step :: *float (Default = 0.05)*
          Step to test the different theta values. Recommended 0.1 or 0.05

        Returns:
        ==============
        best_theta :: *float*
          Best m/n value for the profile. The one that icreases the linearity
          for the whole Chi-Zi profile
        """
        best_r2 = 0
        best_theta = 0
        theta_values = np.arange(0, 1, step)
        zi = self._data[::-1, 2]

        for theta in theta_values:
            ai = self._data[::-1, 4] ** theta
            a0 = a0 ** theta
            chi = [0]
            for n in range(len(ai)):
                if n > 0:
                    dx = self._data[n, 3] - self._data[n - 1, 3]
                    chi.append(chi[n - 1] + (a0 * dx / ai[n]))

            # Regresion into chi-elevation space to get r^2
            a1 = np.array([chi, np.ones(len(chi))]).T
            y1 = zi
            model, resid = np.linalg.lstsq(a1, y1)[:2]
            r2 = 1 - resid / (y1.size * y1.var())
            if r2 > best_r2 and best_theta != 0:
                best_r2 = r2
                best_theta = theta

        return best_theta


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


# Debug
DEM = "../../test/data/darro25.tif"
Flow_accumulation = "../../test/data/darro25fac.tif"
Output_channels = "../../test/data/out_main_ch.shp"
Use_threshold = False
Threshold = 1000
Units = "CELL"
Use_basins = False
Basin_shapefile = "../../test/data/cuencas.shp"
Use_heads = True
Heads_shapefile = "../../test/data/main_channels.shp"
Id_field = "id"
# End debug

main(DEM, Flow_accumulation, Use_threshold, Threshold, Units, Use_basins, Basin_shapefile,
     Use_heads, Heads_shapefile, Id_field, Output_channels)