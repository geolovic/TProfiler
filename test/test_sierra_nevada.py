# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 09 March 2017
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
    fac = "D:/Usuarios/Vicente/Desktop/sierranevada_santamarta_Granada/FlowAcc_projected_weightRainfall/FlowAcc_projected_weight.tif"
    dem = "D:/Usuarios/Vicente/Desktop/sierranevada_santamarta_Granada/DEM_SNSM/DEM_proj_fill_SNSM.tif"
    umbral = 25000
    units = "CELL"

    cabeceras = p.get_heads(fac, dem, umbral, units=units)
    perfiles = p.get_profiles(fac, dem, cabeceras, tributaries=True)

    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 40
    out_txt = "D:/Usuarios/Vicente/Desktop/perfiles.txt"
    outfile = open(out_txt, "w")
    outfile.write("perfil;X;Y;Chi\n")
    id_perfil = 0
    for perfil in perfiles:
        xi = perfil.get_x().astype("str")
        yi = perfil.get_y().astype("str")
        chi = perfil.get_chi()[::-1].astype("str")
        for n in range(len(xi)):
            linea = str(id_perfil) + ";" + xi[n] + ";" + yi[n] + ";" + chi[n] + "\n"
            outfile.write(linea)
        id_perfil += 1
    outfile.close()


def draw_profiles(perfiles):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(20, 6))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Get the longest profile (i.e. the first one)
    max_l = perfiles[0].length()

    for perfil in perfiles:
        chi, zi = perfil.get_z_chi()
        li = perfil.get_l(False)
        li = max_l - li
        xi = perfil.get_x()
        yi = perfil.get_y()
        z = perfil.get_z()
        ax1.plot(xi, yi, color="b")
        ax2.plot(li, z, color="b")
        ax3.plot(chi, zi, color="b")

    perfil = perfiles[0]
    chi, zi = perfil.get_z_chi()
    li = perfil.get_l(False)
    li = max_l - li
    xi = perfil.get_x()
    yi = perfil.get_y()
    z = perfil.get_z()
    ax1.plot(xi, yi, color="r")
    ax2.plot(li, z, color="r")
    ax3.plot(chi, zi, color="r")
    plt.show()


test01()
