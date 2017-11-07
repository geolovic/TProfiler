# -*- coding: iso-8859-15 -*-
#
#  ProfilerApp.py
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

#  Version: 1.2
#  November, 6th 2017

#  Last modified November, 6th 2017

import matplotlib.pyplot as plt
import numpy as np
import ogr
import osr
import os
import pickle


# QGIS TOOLBOX CODE
# ===================
##Input_File=file
##Slope_Regression_points=number 0
##Ksn_Regression_points=number 0


# ARGUMENTS
# ===============
in_file = str(Input_File)
slope_reg_points = int(Slope_Regression_points)
ksn_reg_points = int(Ksn_Regression_points)


# Disable some keymap characters that interfere with graph key events
plt.rcParams["keymap.xscale"] = [""]
plt.rcParams["keymap.yscale"] = [""]
plt.rcParams["keymap.save"] = [u'ctrl+s']
key_dict = {"LEFT": -1, "RIGHT": 1, "+": 1, "-": -1}
k_type = {0: "r*", 1: "g*", 2: "r.", 3: "y."}  # Knickpoint types


class ProfilerApp:
    """
    Class to draw longitudinal, area-slope, chi and ksn profiles.
    Allows selecting knickpoints and do regressions on chi-elevation profiles
    """
    gr_types = {1: 'longitudinal', 2: "area-slope", 3: "chi", 4: "ksn"}

    def __init__(self, figure, profiles, basedir=""):
        """
        :param figure: Matplotlib.Figure to draw graphics
        :param profiles: numpy.array with TProfiles objects
        :param basedir: Base directory to save regressions and knickpoints
        """
        # Create the output folder (for knickpoints and regressions)
        if basedir[-1] == "/":
            basedir = basedir[:-1]
        self.basedir = basedir + "/out_files"
        if not os.path.exists(self.basedir):
            os.mkdir(self.basedir)

        # Variable to avoid double click issue in QGIS for Mac
        self.last_clicked = -1
        
        self.figure = figure
        self.ax = figure.add_subplot(111)
        self.g_type = 1
        self.profiles = profiles
        self.active = 0
        self.mode = None
        self.regression = None
        self.mode_txt = self.figure.text(0.025, 0.95, "", color="navy")
        self.n_profiles = len(profiles)
        self.knick_points = [[] for n in range(self.n_profiles)]
        self.regressions = [[] for n in range(self.n_profiles)]
        self.reg_points = []  # Temporary list of 2 elements [ind0, ind1]
        self.dam_points = []  # Temporary list of 2 elements [ind0, ind1]
        self.activate()
        self.draw()

    def activate(self):
        self.kcid = self.figure.canvas.mpl_connect("key_press_event", self.key_input)

    # DRAWING METHODS
    # ===============
    def draw(self):
        """
        Draw the selected graphic for the active profile
        """
        if self.g_type == 1:
            self._draw_long_profile()
        elif self.g_type == 2:
            self._draw_area_slope()
        elif self.g_type == 3:
            self._draw_chi_profile()
        elif self.g_type == 4:
            self._draw_ksn_profile()
        self._mode_inf()
        self.figure.canvas.draw()

    def _draw_long_profile(self):
        """
        Draw the longitudinal profile for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        li = perfil.get_l() / 1000.
        zi = perfil.get_z()
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(li, zi, c="#1f66b4", lw=1., picker=4)
        
        # Some parameters to make graphic nicer
        dl = (li.max() - li.min()) * 0.05
        dz = (zi.max() - zi.min()) * 0.05
        self.ax.set_xlim(li.min() - dl, li.max() + dl)
        self.ax.set_ylim(zi.min() - dz, zi.max() + dz)

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(li[k[0]], zi[k[0]], k_type[k[1]], mew=0.7, mec="k", ms=12)

    def _draw_area_slope(self):
        """
        Draw area-slope log-log plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        slopes = perfil.get_slope()
        areas = perfil.get_area()
        self.ax.set_xlabel("Area $m^2$")
        self.ax.set_ylabel("Slope (reg. points: {0})".format(perfil.slope_reg_points))
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.set_title(perfil.name)
        self.ax.plot(areas, slopes, "b+", mew=0.5, picker=4)
        
        # Some parameters to make graphic nicer
        da = (np.log(areas.max()) - np.log(areas.min())) * 0.05
        ds = (np.log(slopes.max()) - np.log(slopes.min())) * 0.05
        self.ax.set_xlim(np.exp(np.log(areas.min()) - da), np.exp(np.log(areas.max()) + da))
        self.ax.set_ylim(np.exp(np.log(slopes.min()) - ds), np.exp(np.log(slopes.max()) + ds))
        
        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(areas[k[0]], slopes[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=12)

    def _draw_chi_profile(self):
        """
        Draw chi-elevation plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        chi = perfil.get_chi()
        zi = perfil.get_z()
        self.ax.set_xlabel("Chi [m]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(chi, zi, c="#1f77b4", lw=1., picker=4)
        
        # Some parameters to make graphic nicer
        dchi = (chi.max() - chi.min()) * 0.05
        dz = (zi.max() - zi.min()) * 0.05
        self.ax.set_xlim(chi.min() - dchi, chi.max() + dchi)
        self.ax.set_ylim(zi.min() - dz, zi.max() + dz)

        # Draw regressions
        for r in self.regressions[self.active]:
            self.ax.plot(r[4][:, 0], r[4][:, 1], c="r", ls="--", lw=1.)
            p1x = self.profiles[self.active].get_chi()[r[0]]
            p1y = self.profiles[self.active].get_z()[r[0]]
            p2x = self.profiles[self.active].get_chi()[r[1]]
            p2y = self.profiles[self.active].get_z()[r[1]]
            self.ax.plot([p1x, p2x], [p1y, p2y], "k+",  mew = 1., markersize=8)

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(chi[k[0]], zi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=12)

    def _draw_ksn_profile(self):
        """
        Draw ksn-distance plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        ksn = perfil.get_ksn()
        li = perfil.get_l() / 1000.
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Ksn (reg. points: {0})".format(perfil.ksn_reg_points))
        self.ax.set_title(perfil.name)
        self.ax.plot(li, ksn, c="#1f77b4", lw=1., picker=4)

        # Some parameters to make graphic nicer
        dl = (li.max() - li.min()) * 0.05
        dk = (ksn.max() - ksn.min()) * 0.05
        self.ax.set_xlim(li.min() - dl, li.max() + dl)
        self.ax.set_ylim(ksn.min() - dk, ksn.max() + dk)
        
        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(li[k[0]], ksn[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=12)

    def _mode_inf(self):
        """
        Set the text to inform about the drawing mode
        R - Regresion mode
        K - Knickpoint mode
        D - Dam remover mode
        """
        if self.mode:
            self.mode_txt.set_text("Mode: {0}".format(self.mode))
        else:
            self.mode_txt.set_text("")

    # EVENT METHODS
    # =============
    def key_input(self, event):
        """
        Controls keyboard events to add functionality to application
        :param event: keypress event

        Keyboard actions:
            1-4: Changes the graphic type [long, area/slope, chi, ksn]
            R: Enters in regression mode
            K: Enters in knickpoint mode
            D: Enters in dam remover mode
            C: Clears the last regression or last knickpoint (if in mode R or K)
            + / - : Increases/reduces reg points in area-slope or ksn plots
            S: Save all the knickpoints and regressions
            M: Draws a map with all the profiles and knickpoints (hightlighting the active one)
        """
        key = event.key.upper()
        if key in ["1", "2", "3", "4"]:
            self.g_type = int(key)
            self.draw()

        elif key in ["LEFT", "RIGHT"]:
            self.next_profile(key_dict[key])

        elif key in ["K", "R", "D"]:
            self.change_mode(key)

        elif key == "UP":
            self.change_knickpoint_type()

        elif key == "C":
            self.remove_last()

        elif key in ["+", "-"] and self.g_type in [2, 4]:
            self.change_reg_points(key_dict[key])

        elif key == "S":
            self.save()

        # elif key == "M":
        #     self.draw_map()

    def change_knickpoint_type(self):
        if self.mode == "K" and len(self.knick_points[self.active]) > 0:
            kp = self.knick_points[self.active][-1]
            kp[1] = (kp[1] + 1) % len(k_type)
            # self.knick_points[self.active][-1][1] = (self.knick_points[self.active][-1][1] + 1) % 2
            self.draw()

    def pick_point(self, event):
        """
        Pick a point in the active profile
        :param event: matplotlib picker event
        """
        # In the case that many points are picked (true if the profile has several points). Take the center one.
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]
        
        # Code to avoid multiple points clicked at the same position in QGIS for Mac (magic mouse issue??)
        if ind == self.last_clicked:
            return
        else:
            self.last_clicked = ind
        
        # Select the proper option according with the drawing mode
        if self.mode == "K":
            self.knick_points[self.active].append([ind, 0])
            self.draw()

        elif self.mode == "R" and self.g_type == 3:
            self.reg_points.append(ind)
            chi = self.profiles[self.active].get_chi()[ind]
            zi = self.profiles[self.active].get_z()[ind]
            self.ax.plot(chi,  zi,  "k+",  mew = 1., markersize=8)
            self.figure.canvas.draw()
            self.create_regression()

        elif self.mode == "D" and self.g_type == 3:
            self.dam_points.append(ind)
            chi = self.profiles[self.active].get_chi()
            zi = self.profiles[self.active].get_z()
            c_point = chi[ind]
            z_point = zi[ind]
            self.ax.plot(c_point, z_point,  "k+",  mew = 1., markersize=8)
            self.figure.canvas.draw()
            self.remove_dam()

        elif self.mode == "D" and self.g_type == 1:
            self.dam_points.append(ind)
            li = self.profiles[self.active].get_l() / 1000.
            zi = self.profiles[self.active].get_z()
            c_point = li[ind]
            z_point = zi[ind]
            self.ax.plot(c_point, z_point, "k+",  mew = 1., markersize=8)
            self.figure.canvas.draw()
            self.remove_dam2()

    def remove_dam(self):
        """
        Removes points that correspond with a dam in the longitudinal profile
        It will take points from self.dam_points[] >> list with two points (positions)
        """
        if len(self.dam_points) < 2:
            return
        elif len(self.dam_points) > 2:
            self.dam_points = []
            self.draw()

        p1 = self.dam_points[0]
        p2 = self.dam_points[1]
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp

        perfil = self.profiles[self.active]

        # Create a straigt-line between both points
        chi = perfil.get_chi()
        zi = perfil.get_z()
        arr_inds = list(range(p1, p2))

        pto1 = (chi[p1], zi[p1])
        pto2 = (chi[p2], zi[p2])

        m = (pto2[1] - pto1[1]) / float(pto2[0] - pto1[0])
        b = pto1[1] - (m * pto1[0])

        xi = chi[arr_inds]
        yi = xi * m + b

        perfil._data[p1:p2, 2] = yi
        perfil.calculate_slope(perfil.slope_reg_points)
        perfil.calculate_ksn(perfil.ksn_reg_points)
        self.dam_points = []
        self.draw()

    def remove_dam2(self):
        """
        Removes points that correspond with a dam in the longitudinal profile
        It will take points from self.dam_points[] >> list with two points (positions)
        """
        if len(self.dam_points) < 2:
            return
        elif len(self.dam_points) > 2:
            self.dam_points = []
            self.draw()

        p1 = self.dam_points[0]
        p2 = self.dam_points[1]
        if p1 > p2:
            p1, p2 = p2, p1

        perfil = self.profiles[self.active]

        # Create a straigt-line between both points
        li = perfil.get_l() / 1000.
        zi = perfil.get_z()
        arr_inds = list(range(p1, p2))

        pto1 = (li[p1], zi[p1])
        pto2 = (li[p2], zi[p2])

        m = (pto2[1] - pto1[1]) / float(pto2[0] - pto1[0])
        b = pto1[1] - (m * pto1[0])

        xi = li[arr_inds]
        yi = xi * m + b

        perfil._data[p1:p2, 2] = yi
        perfil.calculate_slope(perfil.slope_reg_points)
        perfil.calculate_ksn(perfil.ksn_reg_points)
        self.dam_points = []
        self.draw()

    def create_regression(self):
        """
        Creates a regression in Chi-Elevation space
        It will take points from self.reg_points[] >> list with two points (positions)
        """
        if len(self.reg_points) < 2:
            return
        elif len(self.reg_points) > 2:
            self.reg_points = []
            self.draw()
            return

        p1 = self.reg_points[0]
        p2 = self.reg_points[1]
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp

        perfil = self.profiles[self.active]

        # Regression in chi plot
        chi = perfil.get_chi()[p1:p2 + 1]
        zi = perfil.get_z()[p1:p2 + 1]
        poli, scr = np.polyfit(chi, zi, deg=1, full=True)[:2]
        r2 = float(1 - scr / (zi.size * zi.var()))
        xc = np.linspace(chi[0], chi[-1], num=5)
        yc = np.polyval(poli, xc)
        arr = np.array((xc, yc)).T
        self.regressions[self.active].append((p1, p2, poli[0], poli[1], arr, r2))

        # Once the regression is created, force redraw and clear reg_point list
        self.reg_points = []
        self.draw()

    def change_mode(self, mode):
        """
        Change the interactive mode:
            K >> Knickpoint selector
            R >> Regression mode (only works in Chi profile, self.gtype == 3)
            D >> Damm remover (only works in Longitudinal profile, self.gtype==1)
        """
        if self.mode == mode:
            self.mode = None
            self.figure.canvas.mpl_disconnect(self.pcid)
        else:
            self.mode = mode
            self.pcid = self.figure.canvas.mpl_connect("pick_event", self.pick_point)
        self.draw()

    def remove_last(self):
        """
        Remove the last knickpoint or the regression from the active profile
        """
        if self.mode == "K" and len(self.knick_points[self.active]) > 0:
            self.knick_points[self.active].pop()
            self.draw()
        elif self.mode == "R" and len(self.regressions[self.active]) > 0:
            self.regressions[self.active].pop()
            self.draw()

    def change_reg_points(self, increment):
        """
        Changes the regressions points in slope-area and ksn profile
        """
        perfil = self.profiles[self.active]
        if self.g_type == 2:
            reg_points = perfil.slope_reg_points + increment
            if reg_points > 2:
                perfil.calculate_slope(reg_points)
        elif self.g_type == 4:
            reg_points = perfil.ksn_reg_points + increment
            if reg_points > 2:
                perfil.calculate_ksn(reg_points)

        self.draw()

    def next_profile(self, idx):
        """
        Get the next profile from profile lists and set it as active profile
        """
        self.active += idx
        self.active %= self.n_profiles
        self.draw()

    def draw_map(self):
        """
        Draws a plot with all the analyzed profiles and highlight the active one
        It also draws knickpoints.
        """
        # Create a new Figure and Axes
        figure = plt.figure()
        ax = figure.add_subplot(111)

        # Draw a map with all the profiles
        for idx, perfil in enumerate(self.profiles):
            xi = perfil.get_x()
            yi = perfil.get_y()
            ax.plot(xi, yi, c="0.6", lw=0.5)
            # Draw the nickpoints
            for k in self.knick_points[idx]:
                ax.plot(xi[k[0]], yi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=7)

        # Draw the active profile with different color
        perfil = self.profiles[self.active]
        xi = perfil.get_x()
        yi = perfil.get_y()
        ax.plot(xi, yi, c="c", lw=1)
        # Draw knickpoints of the active profile bigger
        for k in self.knick_points[self.active]:
            ax.plot(xi[k[0]], yi[k[0]], k_type[k[1]], mew=0.5, mec="k", ms=10)

    # SAVING METHODS
    # =============
    def _save_regressions(self, out_file):

        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.profiles[self.active].get_projection())
        driver = ogr.GetDriverByName("ESRI Shapefile")

        # Check if dataset already exists, create dataset and get layer
        if os.path.exists(out_file):
            dataset = driver.Open(out_file, 1)
            dataset.DeleteLayer(0)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)
        else:
            dataset = driver.CreateDataSource(out_file)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)

        # Create fields
        layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn("r2", ogr.OFTReal))

        # Save the regressions
        for idx, perfil in enumerate(self.profiles):
            for reg in self.regressions[idx]:
                xi = perfil.get_x()[reg[0]:reg[1]]
                yi = perfil.get_y()[reg[0]:reg[1]]
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField("ksn", reg[2])
                feat.SetField("r2", reg[5])

                # Creamos geometria
                geom = ogr.Geometry(ogr.wkbLineString)
                for n in range(len(xi)):
                    geom.AddPoint(xi[n], yi[n])
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

    def _save_knickpoints(self, out_file):
        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.profiles[self.active].get_projection())
        driver = ogr.GetDriverByName("ESRI Shapefile")

        # Check if dataset already exists, create dataset and get layer
        if os.path.exists(out_file):
            dataset = driver.Open(out_file, 1)
            for n in range(dataset.GetLayerCount()):
                dataset.DeleteLayer(n)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)
        else:
            dataset = driver.CreateDataSource(out_file)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)

        # Create fields
        layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))

        # Save knickpoints
        for idx, perfil in enumerate(self.profiles):
            for kp in self.knick_points[idx]:
                xi = perfil.get_x()[kp[0]]
                yi = perfil.get_y()[kp[0]]
                ksn = perfil.get_ksn()[kp[0]]
                tipo = int(kp[1])

                # Populate the field ksn
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField("ksn", ksn)
                feat.SetField("type", tipo)

                # Create geometry
                geom = ogr.Geometry(ogr.wkbPoint)
                geom.AddPoint(xi, yi)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

    def save(self):
        """
        Saves selected knickpoints and regressions into shapefiles
        """
        out_knicks = self.basedir + "/knickpoints.shp"
        out_regres = self.basedir + "/regressions.shp"
        self._save_regressions(out_regres)
        self._save_knickpoints(out_knicks)

        if self.mode:
            self.mode = None
            self.figure.canvas.mpl_disconnect(self.pcid)

        # np.save(self.basedir + "/graph.npy", np.array([self]))
        # np.save(self.basedir + "/graph_profiles.npy", self.profiles)
        self.mode_txt.set_text('Files saved on "{0}"'.format(self.basedir))
        self.ax.figure.canvas.draw()


# IMPORTED MODULES (November, 6th 2017)
# =====================================
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

    def __init__(self, pf_data, dem_res=0, rid=0, thetaref=0.45, chi0=0, reg_points=4, srs="", name="", mouthdist=0,
                 smooth=0):
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
        self._srs = srs  # WKT with the Spatial Reference
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

    def get_area(self, head=True):
        """
        Returns a numpy.array with drainage area values for all vertices

        :param head: boolean - Specifies if areas are returned from head (True) or mouth (False)
        :return: numpy.array wiht area values for all vertices
        """
        areas = np.copy(self._data[:, 4])

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

            if (y.size * y.var()) == 0:
                r2 = 0
            else:
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

            if (sample_z.size * sample_z.var()) == 0:
                r2 = 0
            else:
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


# PROGRAM CODE
# ============
def load_qgis_data(in_file):
    profile_data = pickle.load(open(in_file, "rb"))
    perfiles = []
    for row in profile_data:
        perfiles.append(TProfile(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], 0))
    return perfiles


perfiles = load_qgis_data(in_file)

# Recalculate slope and ksn if desired
if slope_reg_points:
    for perfil in perfiles:
        perfil.calculate_slope(slope_reg_points)

if ksn_reg_points:
    for perfil in perfiles:
        perfil.calculate_ksn(ksn_reg_points)

basedir = os.path.dirname(in_file)
fig = plt.figure()
pgraph = ProfilerApp(fig, perfiles, basedir)

plt.show()