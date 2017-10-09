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
import numpy as np

print("Tests for profiler.profiles_to_shp()")

    
def test01():
    """
    Test for profiles_to_shp() function
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para profiles_to_shp() function")
    print("Output point and line shapefiles")
    print("Test in progress...")
    
    # Test parameters
    profiles_path = "data/in/profiles_basins.npy"
    out_point_shp = "data/out/05_profile_points" + ".shp"
    dl = 250
    out_line_shp = "data/out/05_profile_lines" + ".shp"

    perfiles = np.load(profiles_path)
    p.profiles_to_shp(out_point_shp, perfiles)
    p.profiles_to_shp(out_line_shp, perfiles, dl)

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


test01()
