# -*- coding: utf-8 -*-
"""
José Vicente Pérez
Granada University (Spain)
October, 2017
Testing suite for profiler.py
Last modified: 27 October 2017
"""

import time
import profiler as p
import ogr
import matplotlib.pyplot as plt
import numpy as np
print("Tests for profiler.get_profiles()")

    
def test01():
    """
    Test for get_profiles() function
    Obtiene todos los perfiles en un DEM
        - No cuencas
        - No cabeceras
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 01 para get_profiles() function")
    print("Getting all profiles in a dem")
    print("Test in progress...")
    
    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    umbral = 1000
    units = "CELL"
    out_txt = "data/out/03_get_profiles_test01.txt"

    cabeceras = p.get_heads(fac, dem, umbral, units=units)
    perfiles = p.get_profiles(fac, dem, cabeceras, tributaries=True)
    save_profiles(perfiles, out_txt)
    # draw_profiles(perfiles)  # Desmarcar para pintar los perfiles
    
    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("Result in " + out_txt)
    print("=" * 40)


def test02():
    """
    Test for get_profiles() function
    Testeando una sola cuenca
    Le damos valores de thetaref, reg_points y smooth
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 02 para get_profiles() function")
    print("Testing only the Darro basin")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    basin_shp = "data/in/cuenca_darro.shp"
    main_heads = "data/in/main_heads.shp"
    id_field = "id"
    umbral = 1000
    units = "CELL"
    out_txt = "data/out/03_get_profiles_test02.txt"

    # Get input basin geometry
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(basin_shp)
    layer = dataset.GetLayer(0)
    feat = layer.GetFeature(0)
    basin = feat.GetGeometryRef()

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_heads, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Obtenemos cabeceras dentro de cuenca
    cabeceras = p.heads_inside_basin(heads, basin)

    # Extraemos los perfiles para esas cabeceras
    perfiles = p.get_profiles(fac, dem, cabeceras, basin, tributaries=True, thetaref=0.5, reg_points=10)
    save_profiles(perfiles, out_txt)
    # draw_profiles(perfiles)  # Desmarcar para pintar los perfiles

    fin = time.time()
    print("Test finalizado en " + str(fin - inicio) + " segundos")
    print("=" * 40)


def test03():
    """
    Test for get_profiles() function
    Testeando todas las cuencas de un shapefile
    """
    inicio = time.time()
    print("=" * 40)
    print("Test 03 para get_profiles() function")
    print("Testing all the basins")
    print("Test in progress...")

    # Test parameters
    fac = "data/in/darro25fac.tif"
    dem = "data/in/darro25.tif"
    basin_shp = "data/in/cuencas.shp"
    main_heads = "data/in/main_heads.shp"
    id_field = "id"
    umbral = 1000
    units = "CELL"
    out_txt = "data/out/03_get_profiles_test03.txt"

    # Obtenemos layer con cuencas
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(basin_shp)
    layer = dataset.GetLayer(0)

    # Obtenemos todas las cabeceras del DEM
    heads = p.get_heads(fac, dem, umbral, units)

    # Obtenemos todas las cabeceras de la capa de puntos
    main_heads = p.heads_from_points(dem, main_heads, id_field=id_field)

    # Combinamos las cabeceras de la capa de puntos con las cabeceras del DEM
    heads = np.append(main_heads, heads, axis=0)

    # Creamos lista vacia de perfiles y la vamos llenado con los perfiles de cada cuenca
    perfiles = []
    for feat in layer:
        # Obtenemos geometria
        basin = feat.GetGeometryRef()

        # Obtenemos cabeceras dentro de cuenca
        cabeceras = p.heads_inside_basin(heads, basin)

        # Extraemos los perfiles para esas cabeceras
        perfiles.extend(p.get_profiles(fac, dem, cabeceras, basin, tributaries=True, thetaref=0.5, reg_points=10))

    # Salvamos perfiles
    save_profiles(perfiles, out_txt)
    # draw_profiles(perfiles)  # Desmarcar para pintar los perfiles
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


def save_profiles(profiles, path):
    """
    Saves a list of TProfiles in a text file. For each profile in the list, all the vertexes are recorded.
    Each record, will have the following columns (separated by ";"):
    'id_profile', 'x', 'y', 'length', 'area', 'slope', 'chi', 'ksn'

    :param profiles: list - List of TProfile objects
    :param path: str - Full path to save the TProfile objects
    :return: None
    """
    out_file = open(path, "w")
    out_file.write("id_profile;x;y;length;area;slope;chi;ksn\n")
    id_profile = 0
    for perfil in profiles:
        xi = perfil.get_x()
        yi = perfil.get_y()
        area = perfil.get_area()
        li = perfil.get_l()
        slope = perfil.get_slope()
        chi = perfil.get_chi()
        ksn = perfil.get_ksn()
        for n in range(perfil.n_points):
            data = [str(id_profile), str(xi[n]), str(yi[n]), str(li[n]),
                    str(area[n]), str(slope[n]), str(chi[n]), str(ksn[n])]
            linea = ";".join(data)
            linea += "\n"
            out_file.write(linea)
        id_profile += 1

    out_file.close()


test01()
test02()
test03()
