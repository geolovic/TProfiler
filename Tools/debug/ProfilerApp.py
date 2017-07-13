# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

import matplotlib.pyplot as plt
import numpy as np
import ogr, osr
import os
from profiler import TProfile

# Disable some keymap characters that interfere with graph key events
plt.rcParams["keymap.xscale"] = [""]
plt.rcParams["keymap.yscale"] = [""]
plt.rcParams["keymap.save"] = [u'ctrl+s']

DIRS = {"left": -1, "right": 1}


class ProfilerApp:
    """
    Class to draw longitudinal, chi and area-slope profiles.
    """
    gr_types = {1: 'longitudinal', 2: "area-slope", 3: "chi", 4: "ksn"}
    
    def __init__(self, figure, profiles, basedir=""):
        """
        
        :param figure: Matplotlib.Figure to draw graphics
        :param profiles: numpy.array with TProfiles objects
        :param basedir: Base directory to save regressions and knickpoints
        """
        self.basedir = basedir
        self.figure = figure
        self.ax = fig.add_subplot(111)
        self.g_type = 1
        self.profiles = profiles
        self.active = 0
        self.mode = None
        self.regression = None
        self.mode_txt = self.figure.text(0.9, 0.95, "", color="navy")
        self.n_profiles = len(profiles)
        self.knick_points = [[] for n in range(self.n_profiles)]
        self.regressions = [[] for n in range(self.n_profiles)]
        self.reg_points = []  # List of 3 elements [ind0, ind1, text_pos]
        self.kcid = self.figure.canvas.mpl_connect("key_press_event", self.key_input)
        self.draw()

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
        self.mode_inf()
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
        self.ax.plot(li, zi, c="b", lw=0.7, picker=4)

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(li[k], zi[k], "r*", mew=0.5, mec="k", ms=10)
        
    def _draw_area_slope(self):
        """
        Draw area-slope log-log plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        slopes = perfil.get_slope()
        areas = perfil.get_area(cells=False)
        self.ax.set_xlabel("Area $m^2$")
        self.ax.set_ylabel("Slope (reg. points: {0})".format(perfil.slope_reg_points))
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.set_title(perfil.name)
        self.ax.plot(areas, slopes, "b+", mew=0.5, picker=4)

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(areas[k], slopes[k], "r*", mew=0.5, mec="k", ms=10)

    def _draw_chi_profile(self):
        """
        Draw chi plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        chi = perfil.get_chi()
        zi = perfil.get_z()
        self.ax.set_xlabel("Chi [m]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(chi, zi, c="b", lw=0.7, picker=4)

        # Draw regressions
        for r in self.regressions[self.active]:
            self.ax.plot(r[4][:, 0], r[4][:, 1], c="r", ls="--", lw=1)
            p1x = self.profiles[self.active].get_chi()[r[0]]
            p1y = self.profiles[self.active].get_z()[r[0]]
            p2x = self.profiles[self.active].get_chi()[r[1]]
            p2y = self.profiles[self.active].get_z()[r[1]]
            self.ax.plot([p1x, p2x], [p1y, p2y], "k+")

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(chi[k], zi[k], "r*", mew=0.5, mec="k", ms=10)

    def _draw_ksn_profile(self):
        """
        Draw ksn-elevation plot for the active profile
        """
        self.ax.clear()
        perfil = self.profiles[self.active]
        ksn = perfil.get_ksn()
        li = perfil.get_l() / 1000.
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Ksn (reg. points: {0})".format(perfil.ksn_reg_points))
        self.ax.set_title(perfil.name)
        self.ax.plot(li, ksn, c="b", lw=0.7, picker=4)

        # Draw knickpoints
        for k in self.knick_points[self.active]:
            self.ax.plot(li[k], ksn[k], "r*", mew=0.5, mec="k", ms=10)

    def mode_inf(self):
        """
        Draw a text in the upper-right corner to inform of the drawing mode (R - Regresion mode; K - Knickpoint mode)
        """
        if self.mode:
            self.mode_txt.set_text("Mode: {0}".format(self.mode))
        else:
            self.mode_txt.set_text("")

    def pick_point(self, event):
        """
        Pick a point in the profile and depending on the drawing mode, add a knickpoint or a regression
        :param event: matplotlib picker event
        """
        # In the case that many points are picked (true if the profile has several points). Take the center one.
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]

        # If we are in knickpoint mode, add a knickpoint and force redraw
        if self.mode == "K":
            self.knick_points[self.active].append(ind)
            self.draw()

        # If we are in regression mode, append one point to self.reg_points
        if self.mode == "R" and self.g_type == 3:
            self.reg_points.append(ind)
            chi = self.profiles[self.active].get_chi()[ind]
            zi = self.profiles[self.active].get_z()[ind]
            self.ax.plot(chi, zi, "k+")
            self.figure.canvas.draw()
            self.create_regression()

    def create_regression(self):
        """
        Creates a regression by using the indices in self.reg_points list
        """
        if len(self.reg_points) < 2:
            return

        p1 = self.reg_points[0]
        p2 = self.reg_points[1]
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp
            
        perfil = self.profiles[self.active]
        
        # Regression in chi plot
        chi = perfil.get_chi()[p1:p2+1]
        zi = perfil.get_z()[p1:p2+1]
       
        poli = np.polyfit(chi, zi, deg=1)
        xc = np.linspace(chi[0], chi[-1], num=5)
        yc = np.polyval(poli, xc)
        arr = np.array((xc, yc)).T
        self.regressions[self.active].append((p1, p2, poli[0], poli[1], arr))
        # Once the regression is created, force redraw and clear reg_point list
        self.reg_points = []
        self.draw()
        
    def key_input(self, event):
        """
        Controls keyboard events to add functionality
        :param event: keypress event
        """
        if event.key in ["1", "2", "3", "4"]:
            self.g_type = int(event.key)
            self.draw()

        if event.key.lower() == "s":
            self.save()

        if event.key.lower() == "k":
            if self.mode == "K":
                self.mode = None
                self.figure.canvas.mpl_disconnect(self.pcid)
                self.draw()
            else:
                self.mode = "K"
                self.pcid = self.figure.canvas.mpl_connect("pick_event", self.pick_point)
                self.draw()
        
        if event.key.lower() == "r":
            if self.mode == "R":
                self.mode = None
                self.regression = None
                self.figure.canvas.mpl_disconnect(self.pcid)
                self.draw()
            else:
                self.mode = "R"
                self.regression = None
                self.pcid = self.figure.canvas.mpl_connect("pick_event", self.pick_point)
                self.draw()
        
        if event.key.lower() == "c":
            if self.mode == "K" and len(self.knick_points[self.active]) > 0:
                self.knick_points[self.active].pop()
                self.draw()
            if self.mode == "R" and len(self.regressions[self.active]) > 0:
                self.regressions[self.active].pop()
                self.draw()

        if event.key.lower() == "+" and self.g_type == 2:
            perfil = self.profiles[self.active]
            reg_points = perfil.slope_reg_points + 1
            perfil.calculate_slope(reg_points)
            self.draw()

        if event.key.lower() == "-" and self.g_type == 2:
            perfil = self.profiles[self.active]
            reg_points = perfil.slope_reg_points - 1
            if reg_points < 1:
                return
            perfil.calculate_slope(reg_points)
            self.draw()

        if event.key.lower() == "m":
            self.draw_map()
            
        if event.key in ["left", "right"]:
            self.get_next_profile(DIRS[event.key])
   
    def get_next_profile(self, idx):
        self.active += idx
        self.active %= self.n_profiles
        self.draw()
            
    def draw_map(self):
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
                ax.plot(xi[k], yi[k], "r*", mew=0.5, mec="k", ms=5)    
        
        # Draw the active profile with different color            
        perfil = self.profiles[self.active]
        xi = perfil.get_x()
        yi = perfil.get_y()
        ax.plot(xi, yi, c="c", lw=1)
            
    def save(self):

        out_knicks = self.basedir + "/knickpoints.shp"
        out_regres = self.basedir + "/regressions.shp"

        sp = osr.SpatialReference()
        sp.ImportFromWkt(self.profiles[self.active].get_projection())
        driver = ogr.GetDriverByName("ESRI Shapefile")

        # Save regressions
        if os.path.exists(out_regres):
            dataset = driver.Open(out_regres, 1)
            dataset.DeleteLayer(0)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)
        else:
            dataset = driver.CreateDataSource(out_regres)
            layer = dataset.CreateLayer("regresions", sp, ogr.wkbLineString)

        layer.CreateField(ogr.FieldDefn("ksn", ogr.OFTReal))

        for idx, perfil in enumerate(self.profiles):
            for reg in self.regressions[idx]:
                xi = perfil.get_x()[reg[0]:reg[1]]
                yi = perfil.get_y()[reg[0]:reg[1]]
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField("ksn", reg[2])

                # Creamos geometria
                geom = ogr.Geometry(ogr.wkbLineString)
                for n in range(len(xi)):
                    geom.AddPoint(xi[n], yi[n])
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

        # Save knickpoints
        if os.path.exists(out_knicks):
            dataset = driver.Open(out_knicks, 1)
            for n in range(dataset.GetLayerCount()):
                dataset.DeleteLayer(n)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)
        else:
            dataset = driver.CreateDataSource(out_knicks)
            layer = dataset.CreateLayer("knickpoints", sp, ogr.wkbPoint)

        for idx, perfil in enumerate(self.knick_points):
            for kp in self.knick_points[idx]:
                xi = perfil.get_x()[kp]
                yi = perfil.get_y()[kp]
                feat = ogr.Feature(layer.GetLayerDefn())

                # Creamos geometria
                geom = ogr.Geometry(ogr.wkbLineString)
                geom.AddPoint(xi, yi)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

# ARGUMENTS
# ==========
profiles_file = "../../test/data/river_chi_profiles.npy"
base_dir = "../../test/data"

# CODE
# ==========
perfiles = np.load(profiles_file)
fig = plt.figure()
pgraph = ProfilerApp(fig, perfiles, base_dir)
plt.show()

        
        
        

