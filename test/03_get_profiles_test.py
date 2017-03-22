# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 13 March 2017
"""

import time
import profiler as p

print "Tests for profiler.get_profiles()"

    
def test01():
    """
    Test for get_profiles() function
    """
    inicio = time.time()
    print "=" * 40
    print "Test 01 para get_profiles() function"
    print "Test in progress..."
    
    # Test parameters
    fac = "data/darro25fac.tif"
    dem = "data/darro25.tif"
    basin = "data/cuenca_darro.shp"
    main_ch = "data/darro_main.shp"
    umbral = 1000
    units = "CELL"

    cabeceras = p.heads_inside_basin(fac, dem, basin, umbral, units, main_ch)
    perfiles = p.get_profiles(fac, dem, cabeceras, basin=basin, tributaries=True)
    
    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40
    draw_profiles(perfiles)
    return perfiles


def draw_profiles(perfiles):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (20, 6))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Get the longest profile (i.e. the first one)
    maxL = perfiles[0].Length()

    for perfil in perfiles:
        chi, zi = perfil.get_z_chi()
        li = perfil.get_L(False)
        li = maxL - li
        xi = perfil.get_X()
        yi = perfil.get_Y()
        z = perfil.get_Z()
        ax1.plot(xi, yi, color = "b")
        ax2.plot(li, z, color = "b")
        ax3.plot(chi, zi, color = "b")

    perfil = perfiles[0]
    chi, zi = perfil.get_z_chi()
    li = perfil.get_L(False)
    li = maxL - li
    xi = perfil.get_X()
    yi = perfil.get_Y()
    z = perfil.get_Z()
    ax1.plot(xi, yi, color="r")
    ax2.plot(li, z, color="r")
    ax3.plot(chi, zi, color="r")
    plt.show()
    
p = test01()