# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

import matplotlib.pyplot as plt
import numpy as np
from profiler import TProfile

# Disable some keymap characters that interfere with graph key events
plt.rcParams["keymap.xscale"] = [""]
plt.rcParams["keymap.yscale"] = [""]
plt.rcParams["keymap.save"] = [u'ctrl+s']

DIRS = {"left":-1, "right":1}


class ProfilerApp:
    gr_types = {1: 'longitudinal', 2: "area-slope", 3: "chi", 4: "ksn"}
    
    def __init__(self, figure, profiles):
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
        self.kcid = self.figure.canvas.mpl_connect("key_press_event", self.key_input)
        self.draw()

    def draw(self):
        if self.g_type == 1:
            self._draw_long_profile()
        elif self.g_type == 2:
            self._draw_area_slope()
        elif self.g_type == 3:
            self._draw_chi_profile()
        
        self.mode_inf()   
        self.figure.canvas.draw()

    def _draw_long_profile(self):
        self.ax.clear()
        perfil = self.profiles[self.active]
        li = perfil.get_l() / 1000.
        zi = perfil.get_z()
        self.ax.set_xlabel("Distance [km]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(li, zi, c="b", lw=0.7, picker=4)
        for k in self.knick_points[self.active]:
            self.ax.plot(li[k], zi[k], "r*", mew=0.5, mec="k", ms=10)
        
    def _draw_area_slope(self):
        self.ax.clear()
        perfil = self.profiles[self.active]
        slopes = perfil.get_slope()
        areas = perfil.get_area(cells=False)
        self.ax.set_xlabel("Area $m^2$")
        self.ax.set_ylabel("Slope (reg. points: {0})".format(perfil.slope_reg_points))
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.set_title(perfil.name)
        self.ax.plot(areas, slopes, "b+",mew=0.5, picker=4)
        for k in self.knick_points[self.active]:
            self.ax.plot(areas[k], slopes[k], "r*", mew=0.5, mec="k", ms=10)

    def _draw_chi_profile(self):
        self.ax.clear()
        perfil = self.profiles[self.active]
        chi = perfil.get_chi()
        zi = perfil.get_z()
        self.ax.set_xlabel("Chi [m]")
        self.ax.set_ylabel("Elevation [m]")
        self.ax.set_title(perfil.name)
        self.ax.plot(chi, zi, c="b", lw=0.7, picker=4)
        for k in self.knick_points[self.active]:
            self.ax.plot(chi[k], zi[k], "r*", mew=0.5, mec="k", ms=10)
        
        for r in self.regressions[self.active]:
            self.ax.plot(r[3][:,0], r[3][:,1], c="r", ls="--", lw=1)
            self.ax.text(50, 50, r[2])
    
    def mode_inf(self):
        if self.mode:
            self.mode_txt.set_text("Mode: {0}".format(self.mode))
        else:
            self.mode_txt.set_text("")

    def pick_point(self, event):
        if len(event.ind) > 2:
            ind = (event.ind[-1] + event.ind[0]) // 2
        else:
            ind = event.ind[0]

        if self.mode == "K":
            self.knick_points[self.active].append(ind)
        
        if self.mode == "R":
            if self.regression:
                self.create_regression(self.regression, ind)
                self.regression = None
            else:
                self.regression = ind
        
        self.draw()
    
    def create_regression(self, p1, p2):
        if p1 > p2:
            tmp = p2
            p2 = p1
            p1 = tmp
            
        perfil = self.profiles[self.active]
        
        # Regression in chi plot
        chi = perfil.get_chi()[p1:p2]
        zi = perfil.get_z()[p1:p2]
       
        poli = np.polyfit(chi, zi, deg=1)
        xc = np.linspace(chi[0], chi[-1], num=5)
        yc = np.polyval(poli, xc)
        arr = np.array((xc, yc)).T
        self.regressions[self.active].append((p1, p2, poli[0], arr))                                
        
    def key_input(self, event):

        if event.key in ["1", "2", "3"]:
            self.g_type = int(event.key)
            self.draw()
                
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
                
        if event.key.lower() == "m":
            self.draw_map()
            
        if event.key in ["left", "right"]:
            self.get_next_profile(DIRS[event.key])
   
    def get_next_profile(self, idx):
        self.active += idx
        self.active = self.active % self.n_profiles
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
            

perfiles = np.load(r"C:\Users\Vicente\PycharmProjects\TProfiler\test\data\chi_profiles.npy\chi_profiles.npy")

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
import matplotlib.cm as cmap
cm = cmap.coolwarm

for perfil in perfiles:
   color = cm(perfil.rid/18)
   chi = perfil.get_chi(head=False)
   zi = perfil.get_z(head=False)
   poli = np.polyfit(chi, zi, deg=1)
   ksn = poli[0]
   color = cm(ksn/7)

   ax.plot(chi, zi, color=color, label=perfil.name)

ax.legend()

#pgraph = ProfilerApp(fig, perfiles)

        
        
        

