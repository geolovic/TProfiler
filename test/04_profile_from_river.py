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
    Darro river with tributaries
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
    print("Extracting the Darro Tributaries")
    print("Test in progress...")
    
    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios_darro.shp"
    id_field = "id"
    name_field = "name"

    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, id_field=id_field, name_field=name_field)
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
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios.shp"
    name_field = "name"
    thetaref = 0.35

    # Check with id and name_field
    fig = plt.figure()
    ax = fig.add_subplot(111)
    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, name_field=name_field, thetaref=thetaref)
    for perfil in perfiles:
        chi = perfil.get_chi()
        zi = perfil.get_z(relative=True)
        ax.plot(chi, zi, label=perfil.name)
    ax.legend()
    ax.set_title("Chi profiles with m/n = {0:.2f}".format(thetaref))
    ax.set_xlabel("Chi [m]")
    ax.set_ylabel("Elevation [m]")
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def draw_profiles(perfiles):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for perfil in perfiles:
        chi = perfil.get_chi(False)
        zi = perfil.get_z(False)
        ax.plot(chi, zi, label=perfil.name)

    ax.legend()
    plt.show()


test01()
test02()
