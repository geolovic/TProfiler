# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
April, 2017

Temporal and deleted codes
Last modified: 01 April 2017
"""
import matplotlib.pyplot as plt


def draw_profiles(perfiles):

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


# def get_chi_map(dem, fac, umbral, units="CELL", basin_shp="", main_ch_shp=""):
#     """
#
#     :param dem:
#     :param fac:
#     :param umbral:
#     :param units:
#     :param basin_shp:
#     :param main_ch_shp:
#     :return:
#     """
#     if basin_shp:
#         basin_geometries = []
#         dataset = ogr.Open(basin_shp)
#         layer = dataset.GetLayer(0)
#         for feat in layer:
#             basin_geometries.append(feat.GetGeometryRef())
#
#     if main_ch_shp:
#         main_channel_pts = []
#         dataset = ogr.Open(main_ch_shp)
#         layer = dataset.GetLayer(0)
#         for feat in layer:
#             main_channel_pts.append(feat.GetGeometryRef())
#         sorted_main_channels = []
#         for basin in basin_geometries:
#             for pto in main_channel_pts:
#                 if basin.Contains(pto):
#                     sorted_main_channels.append(pto)
#                     break
