# -*- coding: iso-8859-15 -*-
#
#  profiler.py
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
#  March 02, 2017

#  Last modified 27 March, 2017

import numpy as np
import math
import praster as p
import ogr

PROFILE_DEFAULT = {'name': "",
                   'thetaref': 0.45,
                   'chi0': 0,
                   'smooth_win': 0,
                   'nPoints': 4,
                   'srs': ""}


def get_heads(fac, dem, umbral, units="CELL"):
    """
    Extracts heads from flow accumulation raster based in a threshold

    :param fac: *str* -- Path to the flow acculumation raster
    :param dem: *str* -- Path to the DEM
    :param umbral: *float* -- Threshold for channel initiation (i.e. for head definition)
    :param units: *str -- Units of threshold ("MAP" / "CELL" ) (Default="MAP")
    :return: *list* -- List of tuples (row, col, X, Y, Z, "id")
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

    # Si no hay ningun pixel que sea considerado una cabecera, devolvemos la lista vacia
    if len(heads) == 0:
        output_heads = heads
    else:
        nheads = len(heads)
        # Ordenamos las posiciones por su elevacion
        heads = np.array(heads)
        ind = heads[:, 4].argsort()
        ind = ind[::-1]
        heads = heads[ind]

        # Anadimos columna con ids y devolvemos lista con cabeceras
        ids = np.arange(nheads).astype("float32").reshape(nheads, 1)
        output_heads = np.append(heads, ids, axis=1)

    return output_heads.tolist()


def heads_from_points(dem, point_shp, names_field=""):
    """
    This function creates a list of tuples (row, col, X, Y, "id") from a point shapefile

    :param dem: *str* -- Path to the DEM
    :param point_shp: *str* -- Path to the point shapefile
    :param names_field: *str* -- String with the field with profile names (Default="")
    :return: *list* -- List of tuples (row, col, X, Y, elev, id)
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
        heads.append((cell[0], cell[1], punto[0], punto[1], elev, name))

    return heads


def heads_inside_basin(fac, dem, basin, umbral, units="CELL", main_ch=""):
    # TODO Cambiar está función para que obtenga las cabeceras dentro de la cuenca a partir de cabeceras de entrada
    # heads_inside_basin(heads, basin, main_ch="")
    """
    This function extracts heads that pass a threshold from a flow accumulation raster inside a determined basin

    :param fac: *str* -- Path to the flow accumulation rater
    :param dem: *str* -- Path to the DEM
    :param basin: *str* -- Path to the basin shapefile (it takes the first polygon of the shapefile)
    :param umbral: *float* -- Threshold for channel initiation (i.e. for head definition)
    :param units: *str -- Units of threshold ("MAP" / "CELL" ) (Default="MAP")
    :return: *list* -- List of tuples (row, col, X, Y, "id")
    :param main_ch: *str* -- Path to the point shapefile with the main channel (it takes the first point)
    :return: *list* -- List of tuples (row, col, X, Y, "id")
    """

    # Get all the heads in the DEM
    heads = get_heads(fac, dem, umbral, units)

    if len(heads) == 0:
        return heads

    # Get the basin shapefile and get only heads inside basin shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(basin)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(0)
    basin_polygon = feat.GetGeometryRef()
    basin_heads = []

    # Take only heads that are within basin polygon
    for head in heads:
        pto = ogr.Geometry(ogr.wkbPoint)
        pto.AddPoint(head[2], head[3])
        if basin_polygon.Contains(pto):
            basin_heads.append(head)

    # Ordenamos cabeceras por elevacion
    heads = np.array(basin_heads)
    ind = heads[:, 4].argsort()
    ind = ind[::-1]
    sorted_heads = heads[ind]
    basin_heads = sorted_heads.tolist()

    # Check if a main channel (head) has been defined
    if main_ch:
        demraster = p.open_raster(dem)
        dataset = driver.Open(main_ch)
        layer = dataset.GetLayer(0)
        feat = layer.GetFeature(0)
        geom = feat.GetGeometryRef()
        point = (geom.GetX(), geom.GetY())
        cell = demraster.xy_2_cell(point)
        elev = demraster.get_cell_value(cell)
        basin_heads.insert(0, (cell[0], cell[1], point[0], point[1], elev, 0))

    return basin_heads


def get_profiles(fac, dem, heads, basin="", tributaries=False, **kwargs):
    """
    This function extracts profiles for each head (determined by heads list). It
    completes rivers until the edge of the dem.
    
    Parameters:
    ================
    fac :: *str*
        Path to the flow accumulation raster
    dem :: *str*
        Path to the Digital Elevation Model
    heads :: *tuple*
        Tuple with heads (row, col, X, Y, Z, name)
    tributaries :: *bool*
        Flag to indicate if source distance for tributaries are recorded
        
    kwargs :: *Arguments for profile creation*
        thetaref :: *float* Reference m/n parameter for calculation
        smooth_win :: *float* Window (meters) to smooth each profile with a movil mean
        nPoints :: *int* Number of points down- and upstream for slope calculation in each pixel

    Returns:
    ===============
    profiles :: *list*
        List with output TProfile Objects
    
    """

    # Get facraster and demrasters as pRaster objects
    facraster = p.open_raster(fac)
    demraster = p.open_raster(dem)

    # Get spatial Reference system from dem raster
    srs = demraster.proj

    # Update profile arguments if some of them have been specified
    opt = PROFILE_DEFAULT
    opt.update(kwargs)
    opt['srs'] = srs

    # Creamos praster auxiliar para valores de Chi0 (con valores por defecto de -1)
    chi_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=-1.0)
    # Empty raster to record intial distances
    if tributaries:
        dist_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=0.0)

    # Obtenemos poligono de cuenca de drenaje si se ha especificado
    geom = None
    if basin:
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataset = driver.Open(basin)
        layer = dataset.GetLayer(0)
        feat = layer.GetFeature(0)
        geom = feat.GetGeometryRef()

    # Empty list to record output TProfile objects
    out_profiles = []

    # Process all the rivers (i.e. all the heads)
    for head in heads:

        # List to store tuples with profile data
        profile_data = []

        # Second list to store processed positions (to record Chi values later on)
        positions = []

        chi0 = 0
        dist0 = 0

        # Calculate data for the first point
        pos = (int(head[0]), int(head[1]))
        point = demraster.cell_2_xy(pos)
        area = facraster.get_cell_value(pos) * facraster.cellsize ** 2
        z = demraster.get_cell_value(pos)
        distance = 0

        # Add first point to profile data
        profile_data.append((point[0], point[1], z, distance, area))

        # Add position to position list
        positions.append(pos)

        # 'Mark' Chi raster to indicate that the pixel was processed
        # Chi raster is marked with a value of 0, latter these 0 values will be replaced with Chi values
        chi_raster.set_cell_value(pos, 0)

        # Obtain the next point following the flow direction
        next_pos = facraster.get_flow(pos)

        # Start bucle to obtain vertexes following river's flow
        while next_pos:

            # Data of the new point
            next_point = demraster.cell_2_xy(next_pos)
            z = demraster.get_cell_value(next_pos)
            area = facraster.get_cell_value(next_pos) * facraster.cellsize ** 2
            distance += math.sqrt((next_point[0] - point[0]) ** 2 + (next_point[1] - point[1]) ** 2)

            # Check if point is still inside basin (if specified)
            if geom:
                pto = ogr.Geometry(ogr.wkbPoint)
                pto.AddPoint(next_point[0], next_point[1])
                if not geom.Contains(pto):
                    break

            # Se anaden a las listas
            profile_data.append((next_point[0], next_point[1], z, distance, area))
            positions.append(next_pos)

            # Se comprueba si esta "marcado";
            # si lo esta, coge el valor de chi0 y termina el bucle
            # sino, se marca (con un 0) y se coge el siguiente punto
            if chi_raster.get_cell_value(next_pos) >= 0:
                chi0 = chi_raster.get_cell_value(next_pos)
                # Se coge tambien la distancia
                if tributaries:
                    dist0 = dist_raster.get_cell_value(next_pos)
                break
            else:
                chi_raster.set_cell_value(pos, 0.0)
                point = tuple(next_point)
                pos = tuple(next_pos)
            next_pos = facraster.get_flow(next_pos)

        # Para que el rio sea valido tiene que tener contener al menos tres pixeles
        if len(profile_data) >= 3:
            # Creamos el perfil
            name = head[5]
            perfil = TProfile(np.array(profile_data), facraster.cellsize, name=name, thetaref=opt['thetaref'],
                              chi0=chi0, smooth_win=opt['smooth_win'], npoints=opt['nPoints'], srs=srs, mouthdist=dist0)
            out_profiles.append(perfil)
            # Rellenamos las posiciones procesadas con los valores de chi
            # Y si hemos seleccionado tributaries, tambiÃ©n se marcan las distancias
            n = 0
            distances = perfil.get_l(False)
            for position in positions:
                chi_raster.set_cell_value(position, perfil._data[n, 7])
                if tributaries:
                    dist_raster.set_cell_value(position, distances[n])
                n += 1

    return out_profiles


# def get_chi_map(dem, fac, umbral, units="CELL", basin_shp="", main_ch_shp=""):
#     """
#
#     :param dem:
#     :param fac:
#     :param umbral:
#     :param units:
#     :param basin_shp:
#     :param main_ch_shp:
#     :return:
#     """
#     if basin_shp:
#         basin_geometries = []
#         dataset = ogr.Open(basin_shp)
#         layer = dataset.GetLayer(0)
#         for feat in layer:
#             basin_geometries.append(feat.GetGeometryRef())
#
#     if main_ch_shp:
#         main_channel_pts = []
#         dataset = ogr.Open(main_ch_shp)
#         layer = dataset.GetLayer(0)
#         for feat in layer:
#             main_channel_pts.append(feat.GetGeometryRef())
#         sorted_main_channels = []
#         for basin in basin_geometries:
#             for pto in main_channel_pts:
#                 if basin.Contains(pto):
#                     sorted_main_channels.append(pto)
#                     break

# TODO Crear una funcion que obtenga features de lineas a partir de un perfil
# TODO Crear una funcion que obtenga features de puntos a partir de un perfil

def profiles_to_shp(path, profiles, distance=0):
    # TODO Escribir documentacion para esta funcion
    """

    :param path:
    :param profiles:
    :param distance:
    :return:
    """
    # Creamos un shapefile (puntos o lineas, segun distance)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDatasource(path)
    sp = profiles[0].srs
    if distance > 0:
        layer = dataset.CreateLayer("perfiles", sp, ogr.wkbLineString)
    else:
        layer = dataset.CreateLayer("perfiles", sp, ogr.wkbPoint)

    # Anadimos campos
    layer.CreateField(ogr.FieldDefn("id_profile", ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn("L", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("chi", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("rksn", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("slope", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("rslope", ogr.OFTReal))

    id_perfil = 1
    for profile in profiles:
        xi = profile.get_x()
        yi = profile.get_y()
        li = profiles.get_l()
        ai = profile.get_a()
        slp = profile.get_slope()[0]
        rslp = profile.get_r2()
        chi = profile.get_chi()
        ksn, rksn = profile.get_ksn(n_points=profile.npoints, full=True)
        for n in range(xi.size):
            feat = ogr.Feature(layer.GetLayerDefn())
            geom = ogr.Geometry(ogr.wkbPoint)
            geom.AddPoint(xi[n], yi[n])
            feat.SetGeometry(geom)
            feat.SetField('id_profile', id_perfil)
            feat.SetField('L', li[n])
            feat.SetField('area', ai[n])
            feat.SetField('chi', chi[n])
            feat.SetField('ksn', ksn[n])
            feat.SetField('rksn', rksn[n])
            feat.SetField('slope', slp[n])
            feat.SetField('rslope', rslp[n])
            layer.CreateFeature(feat)
        id_perfil += 1


class TProfile:
    """
    Properties:
    ============================
    self.dem_res :: *float*
      Dem resolution of the Digital elevation model used to retreive area-elevation data
    self.name :: *str*
      Name of the profile
    self.thetaref :: *float*
      Value of m/n used in area-slope calculations
    self._data :: *numpy.array*
      9-column numpy.array with the input data

    =======   ==============================================
    Column    Description
    =======   ==============================================
    c0        X Coordinates of river profile vertex
    c1        Y Coordinates of river profile vertex
    c2        Z Elevation of the river profile vertex
    c3        L Distance to river head
    c4        A Drainage area to the vertex
    c5        slope of river profile in each vertex
    c6        Quality slope, correlation coefficient (r^2) of the slope regression
    c7        Chi (Integral mode)
    c8        Raw Z Elevation of the river profile vertex (used to reset the profile)
    =======   ==============================================
    """

    def __init__(self, pf_data, dem_res=0, name="", thetaref=0.45, chi0=0,
                 smooth_win=0, npoints=4, srs="", mouthdist=0):
        """
        Class that defines a river profile with morphometry capabilities.

        :param pf_data: *numpy array* -- Array with input values
        :param dem_res: *float* -- Resolution of the DEM used to extract profile features
        :param name: *str* -- Name of the profile (label) (Default="")
        :param thetaref: *float* -- Thetaref (m/n) value used to calculate Chi and Ksn indexes (Default=0.45)
        :param chi0: *float* -- Value of chi index for first point (for tributaries) (Default=0)
        :param smooth_win: *float* -- Window size to smooth the river profile (Default=0)
        :param npoints: *int* -- Number of points (at each side) to calculate slope for each vertex (Default=4)
        :param srs: *str* -- Spatial reference system expresed as well knwon text (wkt) (Default="")
        :param mouthdist: *float* -- Distance from profile to the river mouth (for tributaries) (Default=0)

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
        self._index = 0
        self._srs = srs  # EPSG Code of the Spatial Reference
        self._mouthdist = mouthdist
        self.dem_res = float(dem_res)
        self.name = unicode(name)
        self.thetaref = abs(thetaref)

        # Exit if pf_data is empty
        # This creates an empty RiverProfile object
        if len(pf_data) == 0:
            return

            # Get profile data from pf_data array
        slope = np.empty(len(pf_data))
        slope.fill(np.nan)
        r_slope = np.empty(len(pf_data))
        r_slope.fill(np.nan)
        chi = np.empty(len(pf_data))
        chi.fill(np.nan)
        self._data = np.array((pf_data[:, 0],
                               pf_data[:, 1],
                               pf_data[:, 2],
                               pf_data[:, 3],
                               pf_data[:, 4],
                               slope,
                               r_slope,
                               chi,
                               pf_data[:, 2]))

        self._data = self._data.T
        # Create slopes and chi
        # By default, when created, the river profile object is not smoothed. Only
        # peaks and flat segments are fixed.
        self.smooth(smooth_win)
        # Slopes are calculated for a window of 8 pixels (4 on both sides - except for edge values-)
        self.calculate_slopes(npoints)
        # Chi index is created with with an initial chi value (0 if not specified)
        self.create_chi(chi0=chi0)

    def __len__(self):
        return len(self._data)

    def get_projection(self):
        """
        Returns a string with the projection as wkt
        """
        return self._srs

    def set_projection(self, projection):
        """
        Sets the projection of the data

        Parameters:
        ============================
        projection :: *str*
            String with the projection in wkt
        """
        self._srs = projection

    def length(self):
        """
        Returns the total length of the profile
        """
        return self._data[-1, 3]

    def get_x(self):
        """
        Returns a numpy.array with X values for all vertices.
        """
        return np.copy(self._data[:, 0])

    def get_y(self):
        """
        Returns a numpy.array with Y values for all vertices.
        """
        return np.copy(self._data[:, 1])

    def get_z(self):
        """
        Returns a numpy.array with elevation values for all vertices
        If profile hasn't been smoothed, Z values will be the raw elevations
        """
        return np.copy(self._data[:, 2])

    def get_raw_z(self):
        """
        Returns a numpy.array with raw elevation values (not smoothed)
        """
        return np.copy(self._data[:, 8])

    def get_l(self, head=True):
        """
        Returns a numpy.array with distances for all profile vertices

        :param head: *bool* - Define if distances are considered from head (True) or mouth (False).
        If measured from mouth, a initial distance (self._mouthdist) will be added (to account tributaries)
        :return: numpy.array with distances for all vertices (measured from head or mouth)
        """
        river_length = float(self._data[-1, 3])

        if head:
            li = np.copy(self._data[:, 3])
        else:
            li = river_length - self._data[:, 3] + self._mouthdist

        return li

    def get_a(self):
        """
        Returns a numpy.array with drainage area values for all vertices
        """
        return self._data[:, 4]

    def get_slope(self, threshold=0):
        """
        Returns slopes calculated by linear regression

        :param threshold: *float* R^2 threshold. (Slopes with R^2 < threshold will be in lq_slopes array)
        :return: tuple of arrays (slopes, lq_slopes).
         slopes --> numpy.array of slopes with R^2 >= threshold (lq_slopes will receive a np.nan value)
         lq_slopes --> numpy.array of slopes with R^2 < threshold (slopes will receive a np.nan value)
        """
        slopes = []
        lq_slopes = []
        for n in range(len(self._data)):
            if self._data[n, 6] >= threshold:
                slopes.append(self._data[n, 5])
                lq_slopes.append(np.nan)
            else:
                slopes.append(np.nan)
                lq_slopes.append(self._data[n, 5])

        return np.array(slopes), np.array(lq_slopes)

    def get_r2(self):
        """
        Returns a numpy.array with R2 values of slope linear regressions for all vertices
        """
        return np.copy(self._data[:, 6])

    def get_chi(self):
        """
        Returns chi values for all vertices in ascending order.
        """
        return self._data[::-1, 7]

    def get_z_chi(self, relative_z=False):
        """
        This function returns chi-zi values for the profile

        Returns: tuple (chi, zi)
        =============
        chi :: *numpy.Array*
          Chi values for each vertex in ascending order
        zi :: *numpy.Array*
          Elevation values for each vertex
        """
        chi = self._data[::-1, 7]
        zi = self._data[::-1, 2]
        if relative_z:
            min_elev = zi[0]
            zi = zi - min_elev

        return chi, zi

    def smooth(self, window=-1):
        """
        This method smooths the river profile with a movil mean of window size. If window is not specified, it takes
        a smooth windows of 9 cells (4 of each side of the profile vertex)

        To get the raw elevations use get_raw_Z().
        To reset the smoothing use reset().

        Parameters:
        ===================
        window :: *float*
            Window size to smooth the river profile

        """
        if window == -1:
            n_cells = 4
        else:
            n_cells = int(int((window / self.dem_res) + 0.5) / 2)

        for ind in range(len(self._data)):
            low = ind - n_cells
            high = ind + n_cells + 1
            if low < 0:
                low = 0

            elevations = self._data[low:high, 8]
            self._data[ind, 2] = np.nanmean(elevations)

        # Remove peaks and flat segments
        for n in range(len(self._data) - 1):
            if self._data[n + 1, 2] >= self._data[n, 2]:
                self._data[n + 1, 2] = float(self._data[n, 2]) - 0.001

    def reset(self):
        """
        Reset smooth elevations. When reset, smooth elevations will equal to raw elevations
        """
        for n in range(len(self._data)):
            self._data[n, 2] = np.copy(self._data[n, 8])

    def calculate_slopes(self, n_points=4):
        """
        This function calculates slopes for all vertexes by linear regression of distance-elevation data.
        Slopes are stored in column c5 of self._data. Together with slopes, R^2 are calculated (column  c6)

        :param n_points: Number of profile points before and after each vertex to calculate slope.
        :return: None
        """

        for n in range(len(self._data)):
            low = n - n_points
            high = n + n_points

            if low < 0:
                low = 0

            sample_l = self._data[low:high + 1, 3]
            sample_z = self._data[low:high + 1, 2]

            a = np.array([sample_l, np.ones(len(sample_l))]).T
            y = sample_z
            model, resid = np.linalg.lstsq(a, y)[:2]

            r2 = 1 - resid / (y.size * y.var())
            gradient = model[0]

            self._data[n, 6] = abs(r2)

            if abs(gradient) < 0.001:
                self._data[n, 5] = 0.001
            else:
                self._data[n, 5] = abs(gradient)

    def create_chi(self, a0=1, chi0=0.0):
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

        self._data[:, 7] = chi[::-1]

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

    def get_ksn(self, n_points=4, full=False):
        """
        This function calculates ksn for all vertexes by linear regression of chi-elevation data.

        :param n_points: *int* - Number of profile points before and after each vertex to calculate ksn
        :param full: *bool* - Indicates if r^2 array is also returned
        :return: numpy.array with ksn values for all vertexes. If full is true, it returns a tuple of arrays (ksn, r^2)
        """

        ksn = []
        r_2 = []
        for n in range(len(self._data)):
            low = n - n_points
            high = n + n_points

            if low < 0:
                low = 0

            sample_chi = self._data[low:high + 1, 7]
            sample_z = self._data[low:high + 1, 2]

            a = np.array([sample_chi, np.ones(len(sample_chi))]).T
            y = sample_z
            model, resid = np.linalg.lstsq(a, y)[:2]

            r2 = 1 - resid / (y.size * y.var())
            gradient = model[0]

            r_2.append(abs(r2))
            if abs(gradient) < 0.001:
                ksn.append(0.001)
            else:
                ksn.append(abs(gradient))
        if full:
            return np.array(ksn), np.array(r_2)
        else:
            return np.array(ksn)


def version():
    return "Version: 3.0 - 27 March 2017"
