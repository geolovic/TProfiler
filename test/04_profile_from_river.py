# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
June, 2017
Testing suite for profiler.py
Last modified: 18 June 2017
"""

import time
import profiler as p
import ogr
import matplotlib.pyplot as plt
import numpy as np
print("Tests for profiler.profiles_from_rivers()")

    
def test01():
    """
    Test for profiles_from_rivers() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
    print("Extracting the Darro Tributaries")
    print("Test in progress...")
    
    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    river_shapefile = "data/rios_darro.shp"

    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile)
    draw_profiles(perfiles)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test02():
    """
    Test for profiles_from_rivers() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para profiles_from_rivers function")
    print("Extracting main rivers")
    print("Test in progress...")

    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    river_shapefile = "data/rios.shp"

    thetarefs = np.arange(0.15, 0.6, 0.05)
    fig = plt.figure()

    for idx, theta in enumerate(thetarefs):
        perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, thetaref=theta)
        ax = fig.add_subplot(3, 3, idx+1)
        for perfil in perfiles:
            chi = perfil.get_chi()
            zi = perfil.get_z(relative=True)
            ax.plot(chi, zi, label=perfil.rid)
        ax.set_title("n/m = {0:.2f}".format(theta))

    plt.tight_layout()
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def draw_profiles(perfiles):

    fig = plt.figure(figsize=(20, 6))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Get the longest profile (i.e. the first one)
    perfil0 = perfiles.pop(0)
    total_l = perfil0.length()

    for perfil in perfiles:
        chi = perfil.get_chi(False)
        zi = perfil.get_z(False)
        li = total_l - perfil.get_l(False)[::-1]
        z = perfil.get_z(True)
        xi = perfil.get_x()
        yi = perfil.get_y()

        ax1.plot(xi, yi, color="b")
        ax2.plot(li, z, color="b")
        ax3.plot(chi, zi, color="b")

    chi = perfil0.get_chi(False)
    zi = perfil0.get_z(False)
    li = total_l - perfil0.get_l(False)[::-1]
    z = perfil0.get_z(True)
    xi = perfil0.get_x()
    yi = perfil0.get_y()
    ax1.plot(xi, yi, color="r")
    ax2.plot(li, z, color="r")
    ax3.plot(chi, zi, color="r")
    plt.show()


test01()
test02()
