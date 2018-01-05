# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
June, 2017
Testing suite for profiler.py
Last modified: 29 October 2017
"""

import time
import sys
sys.path.append("../src")
import profiler as p
import matplotlib.pyplot as plt
print("Tests for profiler.profiles_from_rivers()")

    
def test01():
    """
    Test for profiles_from_rivers() function
    Testing all rivers in rivers with tributaries
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
    print("Testing all rivers with tributaries")
    print("Test in progress...")
    
    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios.shp"
    id_field = "id_rio"
    name_field = "name"

    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, id_field=id_field, name_field=name_field)
    title = "Test with id and names. tributaries=True"
    draw_profiles(perfiles, title)

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
    print("Testing all rivers with tributaries, no id, no namefield")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios.shp"

    # Check with id and name_field
    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile)
    title = "Test without id and names. tributaries=True"
    draw_profiles(perfiles, title)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test03():
    """
    Test for profiles_from_rivers() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
    print("Testing all rivers without tributaries")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios.shp"
    id_field = "id_rio"
    name_field = "name"

    # Get profiles
    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, id_field, name_field, tributaries=False)
    title = "Test with ids and names. tributaries=False"
    draw_profiles(perfiles, title)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test04():
    """
    Test for profiles_from_rivers() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_from_rivers function")
    print("Testing a river shapefile with disconected rivers")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    river_shapefile = "data/in/rios2.shp"
    id_field = "id_rio"
    name_field = "name"

    # Get profiles
    perfiles = p.profiles_from_rivers(fac, dem, river_shapefile, id_field, name_field, tributaries=True)
    title = "Test shapefile with disconected rivers"
    draw_profiles(perfiles, title)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def draw_profiles(perfiles, title=""):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for perfil in perfiles:
        chi = perfil.get_chi(False)
        zi = perfil.get_z(False)
        ax.plot(chi, zi, label=perfil.name)
    if title:
        ax.set_title(title)
    ax.legend()
    plt.show()


#test01()
#test02()
#test03()
test04()