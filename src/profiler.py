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

#  Last modified 16 June, 2017

import numpy as np
import math
import praster as p
import ogr
import osr


PROFILE_DEFAULT = {'name': "", 'thetaref': 0.45, 'chi0': 0, 'reg_points': 4, 'srs': "", 'smooth': 0}


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

    # Si no hay ningun pixel que sea considerado una cabecera, devolvemos un array vacio
    if len(heads) == 0:
        output_heads = np.array(heads, dtype="float32").reshape((0, 6))
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
        if id_field in fields:
            hid = str(feat[id_field])
        else:
            hid = n + 1
        n += 1
        heads.append([cell[0], cell[1], punto[0], punto[1], elev, hid])

    # Si no hay ningun punto en el shapefile de cabeceras, devolvemos un array vacio
    if len(heads) == 0:
        output_heads = np.array(heads, dtype="float32").reshape((0, 6))
    else:
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

    if len(basin_heads) == 0:
        output_heads =  np.array([], dtype="float32").reshape((0, 6))
    else:
        output_heads = np.array(basin_heads, dtype="float32")
        output_heads[:, 5] = np.arange(output_heads.shape[0]).astype("float32") + float(first)

    return output_heads


def get_profiles(fac, dem, heads, basin=None, tributaries=False, **kwargs):
    """
    This fucntion extracts profiles for all heads in heads list. It will complete river profiles until the edge of
    the dem or until the basin limits (if basin is specified).
    :param fac: str. String with the path to the flow accumulation raster
    :param dem: str. String with the path to the DEM raster
    :param heads: List of tuples (row, col, X, Y, Z, hid)
    :param basin: ogr.Geometry Polygon that represents one basin
    :param tributaries: Boolean Flag to specify if mouth distances will be calculated
    :param kwargs: TProfile arguments for profile creation
    :return: list of TProfile objects
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

    # Auxiliar PRaster object to record values of Chi
    chi_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=-1.0)

    # Auxiliar PRaster object to record values of distance
    if tributaries:
        dist_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=0.0)

    # Empty list to record output TProfile objects
    out_profiles = []

    # If all the heads come from the same basin, all of them will flow into first river
    # there is no need to check if all vertexes are inside the basin, from the second head all will be inside
    first_river = True

    # Process all the heads
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
        area = facraster.get_cell_value(pos)
        if area == 0:
            area = 1
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
            area = facraster.get_cell_value(next_pos)
            distance += math.sqrt((next_point[0] - point[0]) ** 2 + (next_point[1] - point[1]) ** 2)

            # Check if point is still inside basin (if specified and it is not the first head)
            if basin and first_river:
                pto = ogr.Geometry(ogr.wkbPoint)
                pto.AddPoint(next_point[0], next_point[1])
                if not basin.Contains(pto):
                    profile_data.append((next_point[0], next_point[1], z, distance, area))
                    positions.append(next_pos)
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

        # Para que el rio sea valido tiene que tener contener al menos cuatro pixeles
        if len(profile_data) >= 4:
            # Creamos el perfil
            perfil = TProfile(np.array(profile_data), facraster.cellsize, rid=head[5], thetaref=opt['thetaref'],
                              chi0=chi0, reg_points=opt['reg_points'], srs=srs, mouthdist=dist0, smooth=opt['smooth'])
            out_profiles.append(perfil)
            # Rellenamos las posiciones procesadas con los valores de chi
            # Y si hemos seleccionado tributaries, tambiÃ©n se marcan las distancias
            n = 0
            distances = perfil.get_l(False)
            distances = distances[::-1]
            chi = perfil.get_chi(True)
            for position in positions:
                chi_raster.set_cell_value(position, float(chi[n]))
                if tributaries:
                    dist_raster.set_cell_value(position, float(distances[n]))
                n += 1
        first_river = False
    return out_profiles


def profiles_from_rivers(fac, dem, river_shapefile, id_field="", name_field="", tributaries=True, **kwargs):
    """
    This fucntion extracts profiles for all heads in heads list. It will complete river profiles until the edge of
    the dem or until the basin limits (if basin is specified).
    :param fac: str. String with the path to the flow accumulation raster
    :param dem: str. String with the path to the DEM raster
    :param river_shapefile: str. String with the path to the river shapefile
    :param id_field: str. String with the name of the field that contains indexes to order features
    :param name_field: str. String with the name of the field that contains feature's labels
    :param tributaries: boolean. Boolean to set whether tributaries are calculated or not
    :param kwargs: TProfile arguments for profile creation
    :return: list of TProfile objects
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

    # Auxiliar PRaster object to record values of Chi
    chi_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=0.0)

    # Auxiliar PRaster object to record values of distance
    if tributaries:
        dist_raster = p.create_from_template(dem, p.gdal.GDT_Float32, nodata=0.0)

    # Empty list to record output TProfile objects
    out_profiles = []

    # Open river shapefile and get layer
    dataset = ogr.Open(river_shapefile)
    layer = dataset.GetLayer(0)

    # Get features and order them
    layerdef = layer.GetLayerDefn()
    fields = [layerdef.GetFieldDefn(idx).GetName() for idx in range(layerdef.GetFieldCount())]
    n = 0
    feat_list = []
    idx_list = []

    for feat in layer:
        if id_field in fields:
            hid = int(feat[id_field])
        else:
            hid = n
        feat_list.append(feat)
        idx_list.append(hid)
        n += 1

    feat_array = np.array(feat_list)
    idx_array = np.array(idx_list, dtype="int32")
    ind = np.argsort(idx_array)
    feat_array = feat_array[ind]

    id_profile = 1

    # Process all the lines in river shapefile
    for feat in feat_array:

        # Get geometry and points
        geom = feat.GetGeometryRef()
        puntos = geom.GetPoints()

        # List to store tuples with profile data
        profile_data = []

        # Second list to store processed positions (to record Chi values later on)
        positions = []

        # Calculate data for the first point
        point = puntos[0]
        pos = demraster.xy_2_cell(point)
        area = facraster.get_cell_value(pos)
        if area == 0:
            area = 1
        z = demraster.get_cell_value(pos)
        distance = 0

        # Add first point to profile data
        profile_data.append((point[0], point[1], z, distance, area))

        # Add position to position list
        positions.append(pos)

        # Start bucle to obtain next vertexes
        for next_point in puntos[1:]:
            next_pos = demraster.xy_2_cell(next_point)
            z = demraster.get_cell_value(next_pos)
            area = facraster.get_cell_value(next_pos)
            distance += math.sqrt((next_point[0] - point[0]) ** 2 + (next_point[1] - point[1]) ** 2)

            # Se anaden a las listas
            profile_data.append((next_point[0], next_point[1], z, distance, area))
            positions.append(next_pos)

            # Guardamos el punto para siguiente ejecucion del bucle
            point = tuple(next_point)

        next_pos = demraster.xy_2_cell(next_point)
        # Get Chi value for the last point (and distances if tributaries == True)
        chi0 = chi_raster.get_cell_value(next_pos)
        if tributaries:
            dist0 = dist_raster.get_xy_value(next_point)

        if id_field in fields:
            rid = feat[id_field]
        else:
            rid = id_profile

        if name_field in fields:
            name = feat[name_field]
        else:
            name = str(id_profile)

        # Creamos el perfil
        perfil = TProfile(np.array(profile_data), facraster.cellsize, rid=rid, name=name, thetaref=opt['thetaref'],
                          chi0=chi0, reg_points=opt['reg_points'], srs=srs, mouthdist=dist0, smooth=opt['smooth'])
        out_profiles.append(perfil)

        # Rellenamos las posiciones procesadas con los valores de chi
        # Y si hemos seleccionado tributaries, tambien se marcan las distancias
        n = 0
        distances = perfil.get_l(False)
        distances = distances[::-1]
        chi = perfil.get_chi(True)

        for position in positions:
            chi_raster.set_cell_value(position, float(chi[n]))
            dist_raster.set_cell_value(position, float(distances[n]))
            n += 1

        id_profile += 1

    return out_profiles


def profiles_to_shp(path, profiles, distance=0):
    """
    This function save a list of profiles in a shapefile (point or line shapefile, depending of the distance param)

    :param path: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :param distance: float - Distance for the profile segments. If distance = 0, the shapefile will be a point shp
    :return: None
    """
    if distance == 0:
        _profiles_to_points(path, profiles)
    else:
        _profiles_to_lines(path, profiles, distance)


def _profiles_to_points(path, profiles):
    """
    This function save a list of profiles in a point shapefile

    :param path: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :return: None
    """
    # Creates point shapefle
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(profiles[0].get_projection())
    layer = dataset.CreateLayer("perfiles", sp, ogr.wkbPoint)

    # Add fields
    campos = ["id_profile", "L", "area", "z", "chi", "ksn", "rksn", "slope", "rslope"]
    tipos = [0, 2, 12, 2, 2, 2, 2, 2, 2]
    for n in range(len(campos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))

    # Get data from all the profiles
    id_perfil = 1
    for profile in profiles:
        xi = profile.get_x()
        yi = profile.get_y()
        li = profile.get_l()
        ai = profile.get_area()
        zi = profile.get_z()
        chi = profile.get_chi()
        ksn = profile.get_ksn()
        rksn = profile.get_ksn_r2()
        slp = profile.get_slope()
        rslp = profile.get_slope_r2()
        id_arr = np.zeros(xi.size)
        id_arr.fill(id_perfil)
        data = np.array((xi, yi, id_arr, li, ai, zi, chi, ksn, rksn, slp, rslp)).T
        data = data[:-1]  # Remove the last point, because it will have the area of the merging channel

        for row in data:
            row_list = row.tolist()  # To avoid numpy types
            feat = ogr.Feature(layer.GetLayerDefn())
            geom = ogr.Geometry(ogr.wkbPoint)
            geom.AddPoint(row_list[0], row_list[1])
            feat.SetGeometry(geom)

            for idx, element in enumerate(row_list[2:]):
                feat.SetField(campos[idx], element)

            layer.CreateFeature(feat)

        id_perfil += 1


def _profiles_to_lines(path, profiles, distance):
    """
    This function save a list of profiles in a line shapefile

    :param path: str - Full path to the shapefile
    :param profiles: list - List of TProfile objects
    :param distance: float - Distance for the profile segments
    :return: None
    """

    # Creates point shapefle
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    sp.ImportFromWkt(profiles[0].get_projection())
    layer = dataset.CreateLayer("perfiles", sp, ogr.wkbLineString)

    # Add fields
    campos = ["id_profile", "L", "area", "z", "chi", "ksn", "rksn", "slope", "rslope"]
    tipos = [0, 2, 12, 2, 2, 2, 2, 2, 2]
    for n in range(len(campos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))

    id_perfil = 1
    for perfil in profiles:
        # Transformamos la distancia a número de celdas
        # dx es el número de puntos que tendrá cada segmento (número impar)
        cellsize = perfil.dem_res
        if distance == 0:
            dx = len(perfil._data)
        else:
            dx = int((distance / cellsize) + 0.5)
            if dx % 2 == 0:
                dx += 1

        p1 = 0
        p2 = p1 + dx

        # Get profile data from mouth to head
        xi = perfil.get_x(False)
        yi = perfil.get_y(False)
        chi = perfil.get_chi(False)
        zi = perfil.get_z(False)
        li = perfil.get_l()[::-1]
        area = perfil.get_area(False)

        while p1 < len(chi) - 3:
            # Get the data for the segments
            if p2 >= len(chi):
                p2 = len(chi)
            sample_x = xi[p1:p2]
            sample_y = yi[p1:p2]
            sample_chi = chi[p1:p2]
            sample_zi = zi[p1:p2]
            sample_li = li[p1:p2]

            # Ksn by regression
            coef, resid = np.polyfit(sample_chi, sample_zi, deg=1, full=True)[:2]
            ksn = float(abs(coef[0]))
            rksn = float(1 - resid / (sample_zi.size * sample_zi.var()))

            # Slope by regression
            coef, resid = np.polyfit(sample_li, sample_zi, deg=1, full=True)[:2]
            slope = float(abs(coef[0]))
            rslope = float(1 - resid / (sample_zi.size * sample_zi.var()))

            # Get the other values
            mid_pos = p1 + int((p2 - p1) / 2)
            mid_chi = float(chi[mid_pos])
            mid_area = int(area[mid_pos])
            mid_l = float(li[mid_pos])
            mid_z = float(zi[mid_pos])

            # Get geometry for the line
            geom = ogr.Geometry(ogr.wkbLineString)
            for n in range(len(sample_x)):
                geom.AddPoint(float(sample_x[n]), float(sample_y[n]))

            # Create the feature
            feature = ogr.Feature(layer.GetLayerDefn())

            # Add values to feature
            values = [id_perfil, mid_l, mid_area, mid_z, mid_chi, ksn, rksn, slope, rslope]
            for idx, value in enumerate(values):
                feature.SetField(campos[idx], value)

            # Add geometry to feature and feature to layer
            feature.SetGeometry(geom)
            layer.CreateFeature(feature)

            # Last point of the line was p2-1
            p1 = p2 - 1
            p2 = p1 + dx

        id_perfil += 1


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

    def __init__(self, pf_data, dem_res=0, rid=0, thetaref=0.45, chi0=0, reg_points=4, srs="", name="", mouthdist=0, smooth=0):
        """
        Class that defines a river profile with morphometry capabilities.

        :param pf_data: *numpy array* - Array with input values
        :param dem_res: *float* - Resolution of the DEM used to extract profile features
        :param rid: *int* - Profile Identifier
        :param thetaref: *float* - Thetaref (m/n) value used to calculate Chi and Ksn indexes
        :param chi0: *float* - Value of chi index for first point (for tributaries)
        :param name: *str* - Profile name. It will used as the profile label
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
        if name == "":
            self.name = str(rid)
        else:
            self.name = name
        self.thetaref = abs(thetaref)
        self.slope_reg_points = reg_points
        self.ksn_reg_points = reg_points
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
        self.calculate_slope(self.slope_reg_points)
        self.calculate_chi(chi0=chi0)
        self.calculate_ksn(self.ksn_reg_points)

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
        self.slope_reg_points = reg_points
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
        self.ksn_reg_points = reg_points
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


def version():
    return "Version: 4.2 - 19 June 2017"
