# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
March, 2017

Testing suite for profiler.py
Last modified: 19 June 2017
"""

import time
import profiler as p
import praster as pr
import numpy as np
import matplotlib.pyplot as plt
print("Tests for TProfiler methods")

    
def test01():
    """
    Creates a TProfiler from an array with profile_data
    Test for get_x, get_y
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para TProfiler")
    print("Testing functions get_x(), get_y()")
    print("Test in progress...")
    
    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Test 01 get and print x and y arrays
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    xi1 = perfil.get_x(True)
    yi1 = perfil.get_y(True)
    xi2 = perfil.get_x(False)
    yi2 = perfil.get_y(False)
    ax1.plot(xi1, yi1)
    ax2.plot(xi2, yi2)
    ax1.set_title("head = True")
    ax2.set_title("head = False")
    fig.tight_layout()
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test02():
    """
    Creates a TProfiler from an array with profile_data
    Test for get_l, get_z
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para TProfiler")
    print("Testing functions get_l(), get_z()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Test 01 get and print x and y arrays
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    li1 = perfil.get_l(True)
    zi1 = perfil.get_z(True)
    ax1.plot(li1, zi1)
    ax1.set_title("head = True")
    li2 = perfil.get_l(False)
    zi2 = perfil.get_z(False)
    ax2.plot(li2, zi2)
    ax2.set_title("head = False")
    zi3 = perfil.get_z(True, True)
    ax3.plot(li1, zi3)
    ax3.set_title("Relative elevations, head = True")
    zi4 = perfil.get_z(False, True)
    ax4.plot(li2, zi4)
    ax4.set_title("Relative elevations, head = False")
    fig.tight_layout()
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test03():
    """
    Creates a TProfiler from an array with profile_data
    Test for raw_elevations and smooth
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 03 para TProfiler")
    print("Testing functions smooth() and get_raw_z()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Print raw elevations vs peaks removed elevations
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    li = perfil.get_l(True)
    zi = perfil.get_z(True)
    raw_zi = perfil.get_raw_z(True)
    ax1.plot(li, zi, label="Peaks removed")
    ax1.plot(li, raw_zi, label="Raw elevations")
    ax1.set_title("Raw elevations vs peak removed")
    ax1.legend()
    ax1.set_xlim((6850, 8950))
    ax1.set_ylim((950, 1050))

    # Test for smooth function
    distance = 0
    for n in range(5):
        li = perfil.get_l(True)
        zi = perfil.get_z(True)
        perfil.smooth(distance)
        ax2.plot(li, zi, label=str(distance) + " m")
        distance += 50
    ax2.set_title("Smooth with different distances")
    ax2.legend()
    ax2.set_xlim((8000, 9000))
    ax2.set_ylim((950, 1000))
    fig.tight_layout()
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test04():
    """
    Creates a TProfiler from an array with profile_data
    Test for get_area and get_slopes
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 04 para TProfiler")
    print("Testing functions get_area() and get_slopes()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Get slope area and plot in log scale
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim((1000000, 100000000))
        ax.set_ylim((0.001, 1))

    ai = perfil.get_area(True)
    s1 = perfil.get_slope()
    ax1.plot(ai, s1, "b+")
    ax1.set_title("Raw slopes (all)")

    s2 = perfil.get_slope(threshold=0.9)
    ax2.plot(ai, s2, "b+")
    ax2.set_title("Slopes with threshold >= 0.9")

    s3, lq3 = perfil.get_slope(threshold=0.9, lq=True)
    ax3.plot(ai, lq3, "r+")
    ax3.plot(ai, s3, "b+")
    ax3.set_title("Slopes and low quality slopes (threshold 0.9)")

    s4, lq4 = perfil.get_slope(threshold=0.9, lq=True, head=True)
    a2 = perfil.get_area(head=True)
    ax4.plot(a2, lq4, "r+")
    ax4.plot(a2, s4, "b+")

    ax4.set_title("Example 3 with head=True")

    fig.tight_layout(pad=1)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test05():
    """
    Creates a TProfiler from an array with profile_data
    Test for calculate slopes
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 05 para TProfiler")
    print("Testing functions calculate slopes")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    reg_points = 4
    # Get slope area and plot in log scale
    fig = plt.figure(figsize=(12, 6))
    for n in range(1, 9, 2):
        ax1 = fig.add_subplot(4, 2, n)
        ax2 = fig.add_subplot(4, 2, n+1)
        perfil.calculate_slope(reg_points)
        si = perfil.get_slope()
        ai = perfil.get_area()
        ax1.plot(ai, si, "b+")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_xlim((1000000, 100000000))
        ax1.set_ylim((0.001, 1))
        ax1.set_title("reg_points = " + str(reg_points) + " (normal elevations)")

        perfil.calculate_slope(reg_points, True)
        si = perfil.get_slope(0.9)
        ax2.plot(ai, si, "b+")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlim((1000000, 100000000))
        ax2.set_ylim((0.001, 1))
        ax2.set_title("reg_points = " + str(reg_points) + " (raw elevations)")

        reg_points += 4

    fig.tight_layout(pad=1)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test06():
    """
    Creates a TProfiler from an array with profile_data
    Test for calculate_chi() and get_chi()
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 06 para TProfiler")
    print("Testing functions get_chi() and calculate_chi()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Get slope area and plot in log scale
    fig = plt.figure()
    theta = 0.35
    for n in range(1, 10):
        ax = fig.add_subplot(3, 3, n)
        perfil.thetaref = theta
        perfil.calculate_chi()
        chi = perfil.get_chi(False, True)
        zi = perfil.get_z(False, True)
        ax.plot(chi, zi)
        ax.set_title("Thetaref = {0:.2f}".format(theta))
        theta += 0.05

    fig.tight_layout(pad=1)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test07():
    """
    Creates a TProfiler from an array with profile_data
    Test for get_ksn()
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 07 para TProfiler")
    print("Testing function get_ksn()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    # Get slope area and plot in log scale
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    li = perfil.get_l(True)
    ksn1 = perfil.get_ksn()
    ax1.plot(li, ksn1, "b+")
    ax1.set_title("Raw ksn (all)")

    ksn2 = perfil.get_ksn(threshold=0.9)
    ax2.plot(li, ksn2, "b+")
    ax2.set_title("Ksn with threshold >= 0.9")

    ksn3, lq3 = perfil.get_ksn(threshold=0.9, lq=True)
    ax3.plot(li, lq3, "r+")
    ax3.plot(li, ksn3, "b+")
    ax3.set_title("Ksn and low quality ksn (threshold 0.9)")

    ksn4, lq4 = perfil.get_ksn(threshold=0.9, lq=True, head=False)
    l2 = perfil.get_l(head=False)
    ax4.plot(l2, lq4, "r+")
    ax4.plot(l2, ksn4, "b+")

    ax4.set_title("Example 3 with head=False")

    fig.tight_layout(pad=1)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test08():
    """
    Creates a TProfiler from an array with profile_data
    Test for calculate_ksn
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 08 para TProfiler")
    print("Testing functions calculate_ksn()")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    reg_points = 4

    fig = plt.figure(figsize=(12, 6))
    for n in range(1, 9, 2):
        ax1 = fig.add_subplot(4, 2, n)
        ax2 = fig.add_subplot(4, 2, n + 1)
        perfil.calculate_ksn(reg_points)
        ksn = perfil.get_ksn()
        li = perfil.get_l()
        ax1.plot(li, ksn)
        ax1.set_title("KSN with reg_points = " + str(reg_points) + " (normal elevations)")

        perfil.calculate_ksn(reg_points, raw_z=True)
        ksn = perfil.get_ksn()
        ax2.plot(li, ksn)
        ax2.set_title("KSN with reg_points = " + str(reg_points) + " (raw elevations)")

        reg_points += 4

    fig.tight_layout(pad=1)
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test09():
    """
    Creates a TProfiler from an array with profile_data
    Test for calculate_ksn
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 09 para TProfiler")
    print("Testing ksn and SL plots")
    print("Test in progress...")

    # Test parameters
    pf_data = np.load("data/in/darro_pfdata.npy")
    dem = "data/in/darro25.tif"
    demraster = pr.open_raster(dem)
    srs = demraster.proj
    cellsize = demraster.cellsize

    # Creates the profile
    perfil = p.TProfile(pf_data, cellsize, srs=srs)

    reg_points = 12
    fig = plt.figure()
    ax = fig.add_subplot(111)

    perfil.calculate_ksn(reg_points=reg_points)
    perfil.calculate_slope(reg_points=reg_points)
    li = perfil.get_l()
    slope = perfil.get_slope()
    ksn = perfil.get_ksn()
    sl = slope * li

    sl, = ax.plot(li, sl)
    ax.set_ylabel("SL index")
    ax.set_xlabel("Distance (m)")

    twax = ax.twinx()
    ksn, = twax.plot(li, ksn, color="r")
    twax.set_ylabel("Ksn index")
    twax.legend((sl, ksn), ("SL", "ksn"))
    plt.show()

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


test01()
test02()
test03()
test04()
test05()
test06()
test07()
test08()
test09()
