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
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for perfil in perfiles:
        chi = perfil.get_chi()
        zi = perfil.get_z()
        li = perfil.get_l(head=False)
        zi2 = perfil.get_z(head=False)
        ax1.plot(chi, zi)
        ax2.plot(li, zi2)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test02():
    """
    Test for profiles_from_rivers() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
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

test01()
test02()
