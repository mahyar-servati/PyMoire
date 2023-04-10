#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 00:16:45 2022

@author: Mahyar Servati
"""
import numpy as np
import matplotlib.pyplot as plt
import libtetrabz
import PyMoireFunc as pm
import time

#------------------------------------------------ PyMoire
angle_constant=0
# choose one of twist angles and uncomment it
# angle_constant=5 # for 6.01 degree
# angle_constant=6 # for 5.08 degree
# angle_constant=7 # for 4.41 degree
# angle_constant=8 # for 3.89 degree
# angle_constant=9 # for 3.48 degree
# angle_constant=10 # for 3.15 degree
# angle_constant=11 # for 2.87 degree
angle_constant=31 # for 1.05 degree

pi=np.pi                        # pi number             
c=3.25                          # Interlayer distance
# k-path determination [x, y, z, name, number of points]
kinfo = [[0.3333333, 0.6666667, 0.0, "K", 20],
        [0.0000000, 0.0000000, 0.0, "\u0393", 15],
        [0.5000000, 0.5000000, 0.0, "M", 10],
        [0.6666667, 0.3333333, 0.0, "K'", 0]]

# Reference layer lattice constant
A = 2.4599989204
# Reference layer lattice vector
Rlat = np.array([[np.sqrt(3)/2*A, -A/2],
                  [0.0, A]])

# Reference Atomic Position
Ratm = np.array([[0.0, 0.0],
                  [np.sqrt(3)/3*A, 0.0]])

# Reference lattice constant
A = round(np.sqrt(np.dot(Rlat[0],Rlat[0])), 8)


#---------------------------- Commensorate twist angle
angle = pm.comensorate_moire(angle_constant)
teta, MA, M_angle = pm.incomensorate_moire(angle, A)
#---------------------------- Minimum number of atoms for moire super-cell
natm = int((MA/A)*2)
# natm = 20
#---------------------------- Reference Lattice
start_time = time.time()
Ruc = pm.unit_cell(Rlat)
Rrlat = pm.real_to_reciprocal(Rlat)
RBZ = pm.first_brillouin(Rrlat)
RA, RB = pm.ref_latt(Rlat, Ratm, 3, natm)
end_time = time.time()
print('Reference Lattice done in', round(end_time-start_time, 2), "secound")
#---------------------------- Twisted Lattice
start_time = time.time()
Tlat = pm.twist(teta, Rlat)
Tatm = pm.twist(teta, Ratm)
Tuc = pm.unit_cell(Tlat)
Trlat = pm.real_to_reciprocal(Tlat)
TBZ = pm.first_brillouin(Trlat)
# TA = pm.twist(teta, RA)
# TB = pm.twist(teta, RB)
TA, TB = pm.ref_latt(Tlat, Tatm, 3, natm)
end_time = time.time()
print('Twisted Lattice done in', round(end_time-start_time, 2), "secound")
#---------------------------- Moire Lattice
start_time = time.time()
Mlat = pm.moire_vec(MA, M_angle)
Ma3 = np.array([0.0, 0.0, 20.0])
Muc = pm.unit_cell(Mlat)
Mrlat = pm.real_to_reciprocal(Mlat)
MBZ = pm.first_brillouin(Mrlat)
Mkpath_red, hs, hstitle = pm.k_path(kinfo)
end_time = time.time()
print('Moire Lattice done in', round(end_time-start_time, 2), "secound")
#---------------------------- Interlayer distance between reference and twisted atoms
start_time = time.time()


try:
    lA1A2 = pm.open_from_file(f"distances/lA1A2_{round(angle, 2)}.dat", RA)
    lA1B2 = pm.open_from_file(f"distances/lA1B2_{round(angle, 2)}.dat", RA)
    lB1B2 = pm.open_from_file(f"distances/lB1B2_{round(angle, 2)}.dat", RB)
    lB1A2 = pm.open_from_file(f"distances/lB1A2_{round(angle, 2)}.dat", RB)
    lA2A1 = pm.open_from_file(f"distances/lA2A1_{round(angle, 2)}.dat", TA)
    lA2B1 = pm.open_from_file(f"distances/lA2B1_{round(angle, 2)}.dat", TA)
    lB2B1 = pm.open_from_file(f"distances/lB2B1_{round(angle, 2)}.dat", TB)
    lB2A1 = pm.open_from_file(f"distances/lB2A1_{round(angle, 2)}.dat", TB)
except:
    lA1A2 = pm.distance(A, RA, TA)
    lA1B2 = pm.distance(A, RA, TB)
    lB1B2 = pm.distance(A, RB, TB)
    lB1A2 = pm.distance(A, RB, TA)
    lA2A1 = pm.distance(A, TA, RA)
    lA2B1 = pm.distance(A, TA, RB)
    lB2B1 = pm.distance(A, TB, RB)
    lB2A1 = pm.distance(A, TB, RA)
    pm.save_to_file(f"distances/lA1A2_{round(angle, 2)}.dat", lA1A2)
    pm.save_to_file(f"distances/lA1B2_{round(angle, 2)}.dat", lA1B2)
    pm.save_to_file(f"distances/lB1B2_{round(angle, 2)}.dat", lB1B2)
    pm.save_to_file(f"distances/lB1A2_{round(angle, 2)}.dat", lB1A2)
    pm.save_to_file(f"distances/lA2A1_{round(angle, 2)}.dat", lA2A1)
    pm.save_to_file(f"distances/lA2B1_{round(angle, 2)}.dat", lA2B1)
    pm.save_to_file(f"distances/lB2B1_{round(angle, 2)}.dat", lB2B1)
    pm.save_to_file(f"distances/lB2A1_{round(angle, 2)}.dat", lB2A1)

end_time = time.time()
print('Distance calculation done in', round(end_time-start_time, 2), "secound")
#---------------------------- Making a ribbon 
# rw=0.4

# lA1A2 = pm.ribbon(RA, lA1A2, Mlat, rw)
# lA1B2 = pm.ribbon(RA, lA1B2, Mlat, rw)
# lB1B2 = pm.ribbon(RB, lB1B2, Mlat, rw)
# lB1A2 = pm.ribbon(RB, lB1A2, Mlat, rw)

# lA2A1 = pm.ribbon(TA, lA2A1, Mlat, rw)
# lA2B1 = pm.ribbon(TA, lA2B1, Mlat, rw)
# lB2B1 = pm.ribbon(TB, lB2B1, Mlat, rw)
# lB2A1 = pm.ribbon(TB, lB2A1, Mlat, rw)

# RA = pm.ribbon(RA, RA, Mlat, rw)
# RB = pm.ribbon(RB, RB, Mlat, rw)
# TA = pm.ribbon(TA, TA, Mlat, rw)
# TB = pm.ribbon(TB, TB, Mlat, rw)
#---------------------------- Filtering atom positions and distances to Moire unit-cell
start_time = time.time()
# atom positions
uRA = pm.moire_unitcell(RA, RA, Mlat)
uRB = pm.moire_unitcell(RB, RB, Mlat)
uTA = pm.moire_unitcell(TA, TA, Mlat)
uTB = pm.moire_unitcell(TB, TB, Mlat)
# non-twisted to twisted layer
ulA1A2 = pm.moire_unitcell(RA, lA1A2, Mlat)
ulA1B2 = pm.moire_unitcell(RA, lA1B2, Mlat)
ulB1B2 = pm.moire_unitcell(RB, lB1B2, Mlat)
ulB1A2 = pm.moire_unitcell(RB, lB1A2, Mlat)
# # twisted to non-twisted layer
ulA2A1 = pm.moire_unitcell(TA, lA2A1, Mlat)
ulA2B1 = pm.moire_unitcell(TA, lA2B1, Mlat)
ulB2B1 = pm.moire_unitcell(TB, lB2B1, Mlat)
ulB2A1 = pm.moire_unitcell(TB, lB2A1, Mlat)

end_time = time.time()
print('Atoms and distanses are filtered to Moire unitcell in', round(end_time-start_time, 2), "secound")
#---------------------------- Extracting hoppings parametes from DFT/Wannier calculation
start_time = time.time()
# Intra-layer
tA1A1 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/4-4.dat')
tA1B1 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/4-5.dat')
tB1B1 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/5-5.dat')
eA1A1 = tA1A1[0]
hA1A1 = tA1A1[1]
shA1A1 = tA1A1[2]
eB1B1 = tB1B1[0]
hB1B1 = tB1B1[1]
shB1B1 = tB1B1[2]
hA1B1 = tA1B1[0]
shA1B1 = tA1B1[1]
hB1A1 = hA1B1
shB1A1 = shA1B1

tA2A2 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/9-9.dat')
tA2B2 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/9-10.dat')
tB2B2 = pm.intra_hopping('../QE_Displacing_Scrip/AA/out/10-10.dat')
eA2A2 = tA2A2[0]
hA2A2 = tA2A2[1]
shA2A2 = tA2A2[2]
eB2B2 = tB2B2[0]
hB2B2 = tB2B2[1]
shB2B2 = tB2B2[2]
hA2B2 = tA2B2[0]
shA2B2 = tA2B2[1]
hB2A2 = hA2B2
shB2A2 = shA2B2
# Inter-layer
dA1A2, tA1A2 = pm.extract_hopping(0, 1, '../QE_Displacing_Scrip/AA/out/4-9.dat', 21)
dA1B2, tA1B2 = pm.extract_hopping(0, 1, '../QE_Displacing_Scrip/AB/out/4-9.dat', 21)
# Energy maps
dTE, TE = pm.extract_hopping(0, 13.6056980659, '../QE_Displacing_Scrip/AA/out/Total_Energy_AA.dat', 2)
dVdW, VdW = pm.extract_hopping(0, 13.6056980659, '../QE_Displacing_Scrip/AA/out/VdW_Energy_AA.dat', 2)
dFE, FE = pm.extract_hopping(0, 1, '../QE_Displacing_Scrip/AA/out/Fermi_Energy_AA.dat', 2)
dbg, bg = pm.extract_hopping(0, 1, '../QE_Displacing_Scrip/AA/out/HOMO_LUMO_AA.dat', 4)
end_time = time.time()
print('Hopping extraction done in', round(end_time-start_time, 2), "secound")
#---------------------------- Fixing hopping parameter Outliers
start_time = time.time()
try:
    dA1B2, tA1B2 = pm.extract_hopping(0, 1, '../QE_Displacing_Scrip/AB/out/4-9_non_outlier.dat', 2)
except:
    tA1B2 = pm.outlier(dA1B2, tA1B2)
    houtput = np.append(dA1B2[:,:2], tA1B2, axis=1)
    np.savetxt('../QE_Displacing_Scrip/AB/out/4-9_non_outlier.dat', houtput, fmt='%12.6f %12.6f %12.6f')            
tB1A2 = tA1B2
dB1A2 = np.array([[x*-1, y*-1, r] for [x, y, r] in dA1B2])
dB1B2 = dA1A2
tB1B2 = tA1A2 
print('Outliers removed in', round(end_time-start_time, 2), "secound")
#---------------------------- Interpolation and mapping hopping parameters into moire structure
start_time = time.time()

# # Linear Interpolation
# MtA1A2 = pm.matching_interpolate(dA1A2, tA1A2, lA1A2)
# MtA1B2 = pm.matching_interpolate(dA1B2, tA1B2, lA1B2)
# MtB1B2 = pm.matching_interpolate(dA1A2, tA1A2, lB1B2)
# MtB1A2 = pm.matching_interpolate(dA1B2, tA1B2, lB1A2)

# MtA2A1 = pm.matching_interpolate(dA1A2, tA1A2, lA2A1)
# MtA2B1 = pm.matching_interpolate(dA1B2, tA1B2, lA2B1)
# MtB2B1 = pm.matching_interpolate(dA1A2, tA1A2, lB2B1)
# MtB2A1 = pm.matching_interpolate(dA1B2, tA1B2, lB2A1)

# RBF Interpolation
MtA1A2 = pm.matching_Rbf(dA1A2, tA1A2, lA1A2)
MtA1B2 = pm.matching_Rbf(dA1B2, tA1B2, lA1B2)
MtB1B2 = pm.matching_Rbf(dA1A2, tA1A2, lB1B2)
MtB1A2 = pm.matching_Rbf(dB1A2, tB1A2, lB1A2)

MtA2A1 = pm.matching_Rbf(dA1A2, tA1A2, lA2A1)
MtA2B1 = pm.matching_Rbf(dA1B2, tA1B2, lA2B1)
MtB2B1 = pm.matching_Rbf(dA1A2, tA1A2, lB2B1)
MtB2A1 = pm.matching_Rbf(dB1A2, tB1A2, lB2A1)

# # Seting Hoppings to Zero
# MtA1A2 = np.zeros_like(MtA1A2, dtype=float)
# MtA1B2 = np.zeros_like(MtA1B2, dtype=float)
# MtB1B2 = np.zeros_like(MtB1B2, dtype=float)
# MtB1A2 = np.zeros_like(MtB1A2, dtype=float)

# MtA2A1 = np.zeros_like(MtA2A1, dtype=float)
# MtA2B1 = np.zeros_like(MtA2B1, dtype=float)
# MtB2B1 = np.zeros_like(MtB2B1, dtype=float)
# MtB2A1 = np.zeros_like(MtB2A1, dtype=float)


# energies
MTE = pm.energy_interpolate2D(dTE, TE, lA1A2)
MVdW = pm.energy_interpolate2D(dVdW, VdW, lA1A2)
MFE = pm.energy_interpolate2D(dFE, FE, lA1A2)
Mbg = pm.energy_interpolate2D(dbg, bg, lA1A2)

end_time = time.time()
print('Matching hoppings done in', round(end_time-start_time, 2), "secound")
#---------------------------- Filtering hoppings parameters to Moire unit-cell
start_time = time.time()
# non-twisted to twisted layer
uMtA1A2 = pm.moire_unitcell(RA, MtA1A2, Mlat)
uMtA1B2 = pm.moire_unitcell(RA, MtA1B2, Mlat)
uMtB1B2 = pm.moire_unitcell(RB, MtB1B2, Mlat)
uMtB1A2 = pm.moire_unitcell(RB, MtB1A2, Mlat)
# # twisted to non-twisted layer
uMtA2A1 = pm.moire_unitcell(TA, MtA2A1, Mlat)
uMtA2B1 = pm.moire_unitcell(TA, MtA2B1, Mlat)
uMtB2B1 = pm.moire_unitcell(TB, MtB2B1, Mlat)
uMtB2A1 = pm.moire_unitcell(TB, MtB2A1, Mlat)
# energies
uMTE = pm.moire_unitcell(RA, MTE, Mlat)
uMVdW = pm.moire_unitcell(RA, MVdW, Mlat)
uMFE = pm.moire_unitcell(RA, MFE, Mlat)
uMbg = pm.moire_unitcell(RA, Mbg, Mlat)

end_time = time.time()
print('Hoppings are filtered to Moire unitcell in', round(end_time-start_time, 2), "secound")
#---------------------------- Setting intralayer distance
start_time = time.time()
Mlat3D = np.array([[Mlat[0,0], Mlat[0,1],  0.0],
                   [Mlat[1,0], Mlat[1,1],  0.0],
                   [Ma3[0],   Ma3[1], Ma3[2]]])

intra_AA = pm.intra_distance(Rlat, Ratm[0])
intra_BB = intra_AA
intra_AB = pm.intra_distance(Rlat, Ratm[1])
intra_BA = -1*intra_AB

intra_AA_T = pm.intra_distance(Tlat, Tatm[0])
intra_BB_T = intra_AA_T
intra_AB_T = pm.intra_distance(Tlat, Tatm[1])
intra_BA_T = -1*intra_AB_T

intra_AA = np.append(intra_AA,  np.zeros((len(intra_AA),1),dtype=float), axis=1)
intra_AB = np.append(intra_AB,  np.zeros((len(intra_AB),1),dtype=float), axis=1)
intra_BB = np.append(intra_BB,  np.zeros((len(intra_BB),1),dtype=float), axis=1)
intra_BA = np.append(intra_BA,  np.zeros((len(intra_BA),1),dtype=float), axis=1)

intra_AA_T = np.append(intra_AA_T,  np.zeros((len(intra_AA_T),1),dtype=float), axis=1)
intra_AB_T = np.append(intra_AB_T,  np.zeros((len(intra_AB_T),1),dtype=float), axis=1)
intra_BB_T = np.append(intra_BB_T,  np.zeros((len(intra_BB_T),1),dtype=float), axis=1)
intra_BA_T = np.append(intra_BA_T,  np.zeros((len(intra_BA_T),1),dtype=float), axis=1)

#---------------------------- Converting atom positions and distance vectors to 3D vectors
# atom positions
uRA3D = np.append(uRA,  np.zeros((len(uRA),1),dtype=float), axis=1)
uRB3D = np.append(uRB,  np.zeros((len(uRB),1),dtype=float), axis=1)
uTA3D = np.append(uTA, c*np.ones((len(uTA),1),dtype=float), axis=1)
uTB3D = np.append(uTB, c*np.ones((len(uTB),1),dtype=float), axis=1)
# non-twisted to twisted layer
ulA1A2[:,:,2] = c
ulA1B2[:,:,2] = c
ulB1A2[:,:,2] = c
ulB1B2[:,:,2] = c
# # twisted to non-twisted layer
ulA2A1[:,:,2] = -1*c
ulA2B1[:,:,2] = -1*c
ulB2A1[:,:,2] = -1*c
ulB2B1[:,:,2] = -1*c
#---------------------------- Converting atom positions and distance vectors to Reduced coordination
# atom positions
red_uRA = pm.cart_to_red(Mlat3D, np.array([uRA3D]))[0]
red_uRB = pm.cart_to_red(Mlat3D, np.array([uRB3D]))[0]
red_uTA = pm.cart_to_red(Mlat3D, np.array([uTA3D]))[0]
red_uTB = pm.cart_to_red(Mlat3D, np.array([uTB3D]))[0]
# distances
# intralayer
# non-twisted layer
red_intra_AB = pm.cart_to_red(Mlat3D, np.array([intra_AB]))[0]
red_intra_AA = pm.cart_to_red(Mlat3D, np.array([intra_AA]))[0]
red_intra_BA = pm.cart_to_red(Mlat3D, np.array([intra_BA]))[0]
red_intra_BB = pm.cart_to_red(Mlat3D, np.array([intra_BB]))[0]
# twisted layer
red_intra_AB_T = pm.cart_to_red(Mlat3D, np.array([intra_AB_T]))[0]
red_intra_AA_T = pm.cart_to_red(Mlat3D, np.array([intra_AA_T]))[0]
red_intra_BA_T = pm.cart_to_red(Mlat3D, np.array([intra_BA_T]))[0]
red_intra_BB_T = pm.cart_to_red(Mlat3D, np.array([intra_BB_T]))[0]
# interlayer
# non-twisted layer
red_inter_AB = pm.cart_to_red(Mlat3D, ulA1B2)
red_inter_AA = pm.cart_to_red(Mlat3D, ulA1A2)
red_inter_BA = pm.cart_to_red(Mlat3D, ulB1A2)
red_inter_BB = pm.cart_to_red(Mlat3D, ulB1B2)
# # twisted layer
red_inter_AB_T = pm.cart_to_red(Mlat3D, ulA2B1)
red_inter_AA_T = pm.cart_to_red(Mlat3D, ulA2A1)
red_inter_BA_T = pm.cart_to_red(Mlat3D, ulB2A1)
red_inter_BB_T = pm.cart_to_red(Mlat3D, ulB2B1)
#---------------------------- Interlayer Neighbors distance and hoppings extraction
# distances
# non-twisted layer
red_inter_AA, red_inter_sAA, red_inter_ssAA = pm.neigh_extract(red_inter_AA)
red_inter_AB, red_inter_sAB, red_inter_ssAB = pm.neigh_extract(red_inter_AB)
red_inter_BA, red_inter_sBA, red_inter_ssBA = pm.neigh_extract(red_inter_BA)
red_inter_BB, red_inter_sBB, red_inter_ssBB = pm.neigh_extract(red_inter_BB)
# # twisted layer
red_inter_AA_T, red_inter_sAA_T, red_inter_ssAA_T = pm.neigh_extract(red_inter_AA_T)
red_inter_AB_T, red_inter_sAB_T, red_inter_ssAB_T = pm.neigh_extract(red_inter_AB_T)
red_inter_BA_T, red_inter_sBA_T, red_inter_ssBA_T = pm.neigh_extract(red_inter_BA_T)
red_inter_BB_T, red_inter_sBB_T, red_inter_ssBB_T = pm.neigh_extract(red_inter_BB_T)
# hoppings
# non-twisted layer
uMtA1A2, usMtA1A2, ussMtA1A2 = pm.neigh_extract(uMtA1A2)
uMtA1B2, usMtA1B2, ussMtA1B2 = pm.neigh_extract(uMtA1B2)
uMtB1A2, usMtB1A2, ussMtB1A2 = pm.neigh_extract(uMtB1A2)
uMtB1B2, usMtB1B2, ussMtB1B2 = pm.neigh_extract(uMtB1B2)
# # twisted layer
uMtA2A1, usMtA2A1, ussMtA2A1 = pm.neigh_extract(uMtA2A1)
uMtA2B1, usMtA2B1, ussMtA2B1 = pm.neigh_extract(uMtA2B1)
uMtB2A1, usMtB2A1, ussMtB2A1 = pm.neigh_extract(uMtB2A1)
uMtB2B1, usMtB2B1, ussMtB2B1 = pm.neigh_extract(uMtB2B1)
#---------------------------- Indexing coupled orbitals
# number of atoms in each sublattice
norb = [len(uRA), len(uRB), len(uTA), len(uTB)]
# intralayer atoms
# non-twisted layer
hi_intra_AB, ws_intra_AB = pm.hopping_index(red_uRA, red_uRB, red_intra_AB, norb, kind='A1-B1')
hi_intra_AA, ws_intra_AA = pm.hopping_index(red_uRA, red_uRA, red_intra_AA, norb, kind='A1-A1')
hi_intra_BA, ws_intra_BA = pm.hopping_index(red_uRB, red_uRA, red_intra_BA, norb, kind='B1-A1')
hi_intra_BB, ws_intra_BB = pm.hopping_index(red_uRB, red_uRB, red_intra_BB, norb, kind='B1-B1')
# twisted layer
hi_intra_AB_T, ws_intra_AB_T = pm.hopping_index(red_uTA, red_uTB, red_intra_AB_T, norb, kind='A2-B2')
hi_intra_AA_T, ws_intra_AA_T = pm.hopping_index(red_uTA, red_uTA, red_intra_AA_T, norb, kind='A2-A2')
hi_intra_BA_T, ws_intra_BA_T = pm.hopping_index(red_uTB, red_uTA, red_intra_BA_T, norb, kind='B2-A2')
hi_intra_BB_T, ws_intra_BB_T = pm.hopping_index(red_uTB, red_uTB, red_intra_BB_T, norb, kind='B2-B2')
# interlayre atoms
# non-twisted layer
# 1st neighbor
hi_inter_AB, ws_inter_AB = pm.hopping_index(red_uRA, red_uTB, red_inter_AB, norb, kind='A1-B2')
hi_inter_AA, ws_inter_AA = pm.hopping_index(red_uRA, red_uTA, red_inter_AA, norb, kind='A1-A2')         
hi_inter_BA, ws_inter_BA = pm.hopping_index(red_uRB, red_uTA, red_inter_BA, norb, kind='B1-A2')
hi_inter_BB, ws_inter_BB = pm.hopping_index(red_uRB, red_uTB, red_inter_BB, norb, kind='B1-B2')
# 2st neighbor
hi_inter_sAB, ws_inter_sAB = pm.hopping_index(red_uRA, red_uTB, red_inter_sAB, norb, kind='A1-B2')
hi_inter_sAA, ws_inter_sAA = pm.hopping_index(red_uRA, red_uTA, red_inter_sAA, norb, kind='A1-A2')         
hi_inter_sBA, ws_inter_sBA = pm.hopping_index(red_uRB, red_uTA, red_inter_sBA, norb, kind='B1-A2')
hi_inter_sBB, ws_inter_sBB = pm.hopping_index(red_uRB, red_uTB, red_inter_sBB, norb, kind='B1-B2')
# 3st neighbor
hi_inter_ssAB, ws_inter_ssAB = pm.hopping_index(red_uRA, red_uTB, red_inter_ssAB, norb, kind='A1-B2')
hi_inter_ssAA, ws_inter_ssAA = pm.hopping_index(red_uRA, red_uTA, red_inter_ssAA, norb, kind='A1-A2')         
hi_inter_ssBA, ws_inter_ssBA = pm.hopping_index(red_uRB, red_uTA, red_inter_ssBA, norb, kind='B1-A2')
hi_inter_ssBB, ws_inter_ssBB = pm.hopping_index(red_uRB, red_uTB, red_inter_ssBB, norb, kind='B1-B2')
# twisted layer
# # 1st neighbor
hi_inter_AB_T, ws_inter_AB_T = pm.hopping_index(red_uTA, red_uRB, red_inter_AB_T, norb, kind='A2-B1')
hi_inter_AA_T, ws_inter_AA_T = pm.hopping_index(red_uTA, red_uRA, red_inter_AA_T, norb, kind='A2-A1')         
hi_inter_BA_T, ws_inter_BA_T = pm.hopping_index(red_uTB, red_uRA, red_inter_BA_T, norb, kind='B2-A1')
hi_inter_BB_T, ws_inter_BB_T = pm.hopping_index(red_uTB, red_uRB, red_inter_BB_T, norb, kind='B2-B1')
# 2st neighbor
hi_inter_sAB_T, ws_inter_sAB_T = pm.hopping_index(red_uTA, red_uRB, red_inter_sAB_T, norb, kind='A2-B1')
hi_inter_sAA_T, ws_inter_sAA_T = pm.hopping_index(red_uTA, red_uRA, red_inter_sAA_T, norb, kind='A2-A1')         
hi_inter_sBA_T, ws_inter_sBA_T = pm.hopping_index(red_uTB, red_uRA, red_inter_sBA_T, norb, kind='B2-A1')
hi_inter_sBB_T, ws_inter_sBB_T = pm.hopping_index(red_uTB, red_uRB, red_inter_sBB_T, norb, kind='B2-B1')
# 3st neighbor
hi_inter_ssAB_T, ws_inter_ssAB_T = pm.hopping_index(red_uTA, red_uRB, red_inter_ssAB_T, norb, kind='A2-B1')
hi_inter_ssAA_T, ws_inter_ssAA_T = pm.hopping_index(red_uTA, red_uRA, red_inter_ssAA_T, norb, kind='A2-A1')         
hi_inter_ssBA_T, ws_inter_ssBA_T = pm.hopping_index(red_uTB, red_uRA, red_inter_ssBA_T, norb, kind='B2-A1')
hi_inter_ssBB_T, ws_inter_ssBB_T = pm.hopping_index(red_uTB, red_uRB, red_inter_ssBB_T, norb, kind='B2-B1')

end_time = time.time()
print('Hoppings index calculation done in', round(end_time-start_time, 2), "secound")
#---------------------------- Stacking orbital indexes, hoppings and distance vectors
start_time = time.time()
onsite_A1 = pm.onsite(red_uRA, eA1A1, norb, kind='A1')
onsite_B1 = pm.onsite(red_uRB, eB1B1, norb, kind='B1')
onsite_A2 = pm.onsite(red_uTA, eA2A2, norb, kind='A2')
onsite_B2 = pm.onsite(red_uTB, eB2B2, norb, kind='B2')
# intralayer atoms
# non-twisted layer
hoppings_intra_AA = pm.stack(ws_intra_AA, hi_intra_AA, hA1A1, red_intra_AA)
hoppings_intra_AB = pm.stack(ws_intra_AB, hi_intra_AB, hA1B1, red_intra_AB)
hoppings_intra_BA = pm.stack(ws_intra_BA, hi_intra_BA, hB1A1, red_intra_BA)
hoppings_intra_BB = pm.stack(ws_intra_BB, hi_intra_BB, hB1B1, red_intra_BB)
# twisted layer
hoppings_intra_AA_T = pm.stack(ws_intra_AA_T, hi_intra_AA_T, hA1A1, red_intra_AA_T)
hoppings_intra_AB_T = pm.stack(ws_intra_AB_T, hi_intra_AB_T, hA1B1, red_intra_AB_T)
hoppings_intra_BA_T = pm.stack(ws_intra_BA_T, hi_intra_BA_T, hB1A1, red_intra_BA_T)
hoppings_intra_BB_T = pm.stack(ws_intra_BB_T, hi_intra_BB_T, hB1B1, red_intra_BB_T)
# interlayer atoms
# non-twisted layer
# 1st neighbor
hoppings_inter_AA = pm.stack(ws_inter_AA, hi_inter_AA, uMtA1A2, red_inter_AA)
hoppings_inter_AB = pm.stack(ws_inter_AB, hi_inter_AB, uMtA1B2, red_inter_AB)
hoppings_inter_BA = pm.stack(ws_inter_BA, hi_inter_BA, uMtB1A2, red_inter_BA)
hoppings_inter_BB = pm.stack(ws_inter_BB, hi_inter_BB, uMtB1B2, red_inter_BB)
# 2st neighbor
hoppings_inter_sAA = pm.stack(ws_inter_sAA, hi_inter_sAA, usMtA1A2, red_inter_sAA)
hoppings_inter_sAB = pm.stack(ws_inter_sAB, hi_inter_sAB, usMtA1B2, red_inter_sAB)
hoppings_inter_sBA = pm.stack(ws_inter_sBA, hi_inter_sBA, usMtB1A2, red_inter_sBA)
hoppings_inter_sBB = pm.stack(ws_inter_sBB, hi_inter_sBB, usMtB1B2, red_inter_sBB)
# 3st neighbor
hoppings_inter_ssAA = pm.stack(ws_inter_ssAA, hi_inter_ssAA, ussMtA1A2, red_inter_ssAA)
hoppings_inter_ssAB = pm.stack(ws_inter_ssAB, hi_inter_ssAB, ussMtA1B2, red_inter_ssAB)
hoppings_inter_ssBA = pm.stack(ws_inter_ssBA, hi_inter_ssBA, ussMtB1A2, red_inter_ssBA)
hoppings_inter_ssBB = pm.stack(ws_inter_ssBB, hi_inter_ssBB, ussMtB1B2, red_inter_ssBB)
# twisted layer
# # 1st neighbor
hoppings_inter_AA_T = pm.stack(ws_inter_AA_T, hi_inter_AA_T, uMtA2A1, red_inter_AA_T)
hoppings_inter_AB_T = pm.stack(ws_inter_AB_T, hi_inter_AB_T, uMtA2B1, red_inter_AB_T)
hoppings_inter_BA_T = pm.stack(ws_inter_BA_T, hi_inter_BA_T, uMtB2A1, red_inter_BA_T)
hoppings_inter_BB_T = pm.stack(ws_inter_BB_T, hi_inter_BB_T, uMtB2B1, red_inter_BB_T)
# 2st neighbor
hoppings_inter_sAA_T = pm.stack(ws_inter_sAA_T, hi_inter_sAA_T, usMtA2A1, red_inter_sAA_T)
hoppings_inter_sAB_T = pm.stack(ws_inter_sAB_T, hi_inter_sAB_T, usMtA2B1, red_inter_sAB_T)
hoppings_inter_sBA_T = pm.stack(ws_inter_sBA_T, hi_inter_sBA_T, usMtB2A1, red_inter_sBA_T)
hoppings_inter_sBB_T = pm.stack(ws_inter_sBB_T, hi_inter_sBB_T, usMtB2B1, red_inter_sBB_T)
# 3st neighbor
hoppings_inter_ssAA_T = pm.stack(ws_inter_ssAA_T, hi_inter_ssAA_T, ussMtA2A1, red_inter_ssAA_T)
hoppings_inter_ssAB_T = pm.stack(ws_inter_ssAB_T, hi_inter_ssAB_T, ussMtA2B1, red_inter_ssAB_T)
hoppings_inter_ssBA_T = pm.stack(ws_inter_ssBA_T, hi_inter_ssBA_T, ussMtB2A1, red_inter_ssBA_T)
hoppings_inter_ssBB_T = pm.stack(ws_inter_ssBB_T, hi_inter_ssBB_T, ussMtB2B1, red_inter_ssBB_T)
# stacking all and constructing Hamiltonian
hoppings = np.row_stack((
                        onsite_A1,
                        onsite_B1,
                        onsite_A2,
                        onsite_B2,
                        hoppings_intra_AA,
                        hoppings_intra_AB,
                        hoppings_intra_BA,
                        hoppings_intra_BB,
                        hoppings_intra_AA_T,
                        hoppings_intra_AB_T,
                        hoppings_intra_BA_T,
                        hoppings_intra_BB_T,
                        hoppings_inter_AA,
                        hoppings_inter_AB,
                        hoppings_inter_BA,
                        hoppings_inter_BB,
                        hoppings_inter_sAB,
                        hoppings_inter_sAA,
                        hoppings_inter_sBA,
                        hoppings_inter_sBB,
                        hoppings_inter_ssAB,
                        hoppings_inter_ssAA,
                        hoppings_inter_ssBA,
                        hoppings_inter_ssBB,
                        hoppings_inter_AA_T,
                        hoppings_inter_AB_T,
                        hoppings_inter_BA_T,
                        hoppings_inter_BB_T,
                        hoppings_inter_sAB_T,
                        hoppings_inter_sAA_T,
                        hoppings_inter_sBA_T,
                        hoppings_inter_sBB_T,
                        hoppings_inter_ssAB_T,
                        hoppings_inter_ssAA_T,
                        hoppings_inter_ssBA_T,
                        hoppings_inter_ssBB_T
                        ))

end_time = time.time()
print('Stacking hoppings done in', round(end_time-start_time, 2), "secound")
#---------------------------- Solving Hamiltonian to obtain Band-Structure
start_time = time.time()
print('Solving BandStructure Hamiltonian ....')
# Hamiltonian dimansion
dim = sum(norb)
# converting k-path to 3D reduced coordination vectors
# red_Mkpath = pm.cart_to_red_2D(Mrlat, Mkpath)
# red_Mkpath3D = np.append(red_Mkpath, np.zeros((len(Mkpath),1),dtype=float), axis=1)
red_Mkpath3D = Mkpath_red
nkp=len(red_Mkpath3D)
ret_eval=np.zeros((dim,nkp),dtype=float)
# Solving Hamiltonian in k-path
for n,k in enumerate(red_Mkpath3D):
    kpnt=np.array(k)
    ham=np.zeros((dim,dim),dtype=complex)
    # for i in range(dim):
    #     ham[i,i]=eA1A1
    for hopping in hoppings:
        amp=complex(hopping[5])
        i=hopping[3]
        j=hopping[4]
        length=hopping[6]
        phase=-1*np.exp((2.0j)*np.pi*np.dot(kpnt,length))
        amp=amp*phase
        ham[i,j]+=amp
    ham_use=ham
    eval=np.linalg.eigvalsh(ham_use)
    eval=np.array(eval.real,dtype=float)
    args=eval.argsort()
    eval=eval[args]
    eval = np.array(eval,dtype=float)
    ret_eval[:,n]=eval[:]
    pm.progress_bar(n,len(red_Mkpath3D))
int_evals = ret_eval

# Fermi level adjusting
Ef = {6.01:-1.71, 5.09:-1.71, 4.41:-1.71, 3.89:-1.71,
      3.48:-1.71, 3.15:-1.705, 2.88:-1.707, 1.05:-1.717} 
E_Ef = int_evals - Ef[round(angle, 2)]

# Writing Band-structure in file
band_file = open(
    f"Moire_band_{round(angle, 2)}-degree.dat", "w")
for i in range(int_evals.shape[0]):
    for j in range(int_evals.shape[1]):
        band_file.write(str(j) + " " + str(round(int_evals[i,j], 7)) + "\n")
    band_file.write("\n")
band_file.close()
    
end_time = time.time()
print('BandStructure calculation done in', round(end_time-start_time, 2), "secound")
#---------------------------- DOS_Tetrahedron calculation (please comment if no need to calculate)
start_time = time.time()
print('Solving DOS Hamiltonian ....')
# # moire lattice 3D vector
Mrlat3D = np.append(Mrlat,  np.zeros((len(Mrlat),1),dtype=float), axis=1)
bvec = np.append(Mrlat3D, np.array([[0.0,0.0,20.0]]), axis=0)
# # DOS k-mesh density
ng0 = 10
ng = np.array([ng0, ng0, 1])
# # DOS Hamiltonian
eig = np.empty([ng[0], ng[1], ng[2], dim], dtype=np.float_)
# # Solving Hamiltonian
for i0 in range(ng[0]):
    for i1 in range(ng[1]):
        for i2 in range(ng[2]):
            kvec = np.array([i0, i1, i2]) / ng[0:3]
            ham_k = np.zeros([dim, dim], dtype=np.complex_)
            # for i in range(dim):
            #     ham_k[i,i]=eA1A1
            for hopping in hoppings:
                amp=complex(hopping[5])
                i=hopping[3]
                j=hopping[4]
                length=hopping[6]
                phase=np.exp((2.0j)*np.pi*np.dot(kvec,length))
                amp=amp*phase
                ham_k[i,j]+=amp
            eig[i0, i1, i2, :], ham_k = np.linalg.eigh(ham_k)
    pm.progress_bar(i0,ng[0])
# # Writing eig in file
with open(f"eig_{round(angle, 2)}-degree.dat", "w") as f:
    for i0 in range(eig.shape[0]):
        for i1 in range(eig.shape[1]):
            for i3 in eig[i0,i1,0]:
                f.write(str(round(i3,8)) + " ")
            f.write("\n")

# # Opening saved eig file 
# with open(f"eig_{round(angle, 2)}-degree.dat", "r") as f:
#     lines= f.readlines()
# l=0
# for i0 in range(eig.shape[0]):
#     for i1 in range(eig.shape[1]):
#         eig[i0,i1,0] = np.array(list(float(i) for i in lines[l].split()))
#         l +=1

# # Calculating DOS
eng_div = 600 # number of energy divisions (Its important to DOS accuracy)
energy = np.linspace(-1+Ef[round(angle, 2)], 1+Ef[round(angle, 2)], eng_div)
wght = libtetrabz.dos(bvec, eig, energy)
dos = wght.sum(3).sum(2).sum(1).sum(0)

# # Writing DOS in the file
dos_file = open(f"Moire_dos_{round(angle, 2)}-degree.dat", "w")
for i in range(len(energy)):
    dos_file.write(str(round(energy[i],7)) + " " + str(round(dos[i], 10)) + "\n")
dos_file.close()

# # Calculating fermi level
# ef, wght, iteration = libtetrabz.fermieng(bvec, eig, int(dim/2)-2)
# totalE = wght.sum()
# # writing fermi level to the file
# with open("FE_TE.txt", "a") as f:
#     print(round(angle, 2), round(ef, 8), round(totalE, 8), file=f)

# end_time = time.time()
# print('DOS and Fermi Energy calculation done in', round(end_time-start_time, 2), "secound")

#---------------------------- writing the hr.dat file
# sorting in order to ws and orbital index
hr = hoppings
hr[:,6] = 0
hr = hr[hr[:,3].argsort()]
hr = hr[hr[:,4].argsort(kind='mergesort')]
hr = hr[hr[:,1].argsort(kind='mergesort')]
hr = hr[hr[:,0].argsort(kind='mergesort')]
# ws number
wsnum = len(list(set([(x[0], x[1], x[2]) for x in hr])))
hr = [[x[0],x[1],x[2],x[3]+1,x[4]+1,x[5],x[6]] for x in hr]
# writing to the file
file = open(f'TBG-{round(angle, 2)}-degree_hr.dat', "w")
file.write(f"Tight binding Hamiltonian for Twisted bilayer graphene in twist angle = {round(angle, 2)} degree, written on {time.ctime()}\n")
file.write("         " + str(dim) + "\n")
file.write("         " + str(wsnum) + "\n")
np.savetxt(file, np.ones((wsnum,1)), fmt='%5i', newline='')
file.write("\n")
np.savetxt(file, hr, fmt='%5i %4i %4i %6i %6i %12.6f %4i')
file.close()  
#------------------------------------------------------------- Plotting results
#---------------------------- plot1: Moire unit-cell and BZ
fs=1
fig1=plt.figure(figsize=(10, 10))
fig1.subplots_adjust(wspace=0.1,hspace=0.0)
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = 20

ax1 = fig1.add_subplot(311)
ax1.set_aspect('equal', adjustable='box')
ax1.plot(Ruc[:,0],Ruc[:,1], color='b',linewidth=1)
ax1.plot(Tuc[:,0],Tuc[:,1], color='r',linewidth=1)
ax1.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1.5)
ax1.scatter(RA[:,0], RA[:,1], c='blue', marker=None, s=fs, alpha=0.3) 
ax1.scatter(RB[:,0], RB[:,1], c='blue', marker=None, s=fs, alpha=0.3) 
ax1.scatter(TA[:, 0], TA[:, 1], c='red', marker=None, s=fs, alpha=0.3)
ax1.scatter(TB[:, 0], TB[:, 1], c='red', marker=None, s=fs, alpha=0.3)
ax1.scatter(uRA[:,0], uRA[:,1], c='blue', marker=None, s=fs) 
ax1.scatter(uRB[:,0], uRB[:,1], c='blue', marker=None, s=fs) 
ax1.scatter(uTA[:, 0], uTA[:, 1], c='red', marker=None, s=fs)
ax1.scatter(uTB[:, 0], uTB[:, 1], c='red', marker=None, s=fs)
ax1.arrow(0, 0, Rlat[0,0], Rlat[0,1], head_width = 0.5, width = 0.02, color='b', length_includes_head=True)
ax1.arrow(0, 0, Rlat[1,0], Rlat[1,1], head_width = 0.5, width = 0.02, color='b', length_includes_head=True)
ax1.arrow(0, 0, Tlat[0,0], Tlat[0,1], head_width = 0.5, width = 0.02, color='r', length_includes_head=True)
ax1.arrow(0, 0, Tlat[1,0], Tlat[1,1], head_width = 0.5, width = 0.02, color='r', length_includes_head=True)
ax1.arrow(0, 0, Mlat[0,0], Mlat[0,1], head_width = 0.8, width = 0.02, color='k', length_includes_head=True)
ax1.arrow(0, 0, Mlat[1,0], Mlat[1,1], head_width = 0.8, width = 0.02, color='k', length_includes_head=True)
ax1.set_xlim(xmin=-4.0,xmax=Muc[2,0]+4.0)
ax1.set_ylim(ymin=-4.0,ymax=Muc[2,1]+4.0)
ax1.set_xticks([])
ax1.set_yticks([])

ax2 = fig1.add_subplot(312)
ax2.set_aspect('equal', adjustable='box')
ax2.plot(RBZ[:, 0], RBZ[:, 1], color='b',linewidth=2)
ax2.plot(TBZ[:, 0], TBZ[:, 1], color='r',linewidth=2)
ax2.arrow(0, 0, Rrlat[0,0], Rrlat[0,1], head_width = 0.15, width = 0.02, color='b', length_includes_head=True)
ax2.arrow(0, 0, Rrlat[1,0], Rrlat[1,1], head_width = 0.15, width = 0.02, color='b', length_includes_head=True)
ax2.arrow(0, 0, Trlat[0,0], Trlat[0,1], head_width = 0.15, width = 0.02, color='r', length_includes_head=True)
ax2.arrow(0, 0, Trlat[1,0], Trlat[1,1], head_width = 0.15, width = 0.02, color='r', length_includes_head=True)
LMBZ = np.array(pm.Lattice(10))
for n in LMBZ:
    t = n[0]*Mrlat[0] + n[1]*Mrlat[1]
    if np.sqrt(np.dot(t,t)) < ((np.sqrt(np.dot(Rrlat[0],Rrlat[0])))*((1+0.05)/2)):
        ax2.plot(t[0]+MBZ[:,0], t[1]+MBZ[:,1], color='k',linewidth=1)
ax2.arrow(0, 0, Mrlat[0,0], Mrlat[0,1], head_width = 0.1, width = 0.02, color='k', length_includes_head=True)
ax2.arrow(0, 0, Mrlat[1,0], Mrlat[1,1], head_width = 0.1, width = 0.02, color='k', length_includes_head=True)
ax2.axis('off')


ax3 = fig1.add_subplot(313)
ax3.set_aspect('equal', adjustable='box')
ax3.plot(MBZ[:, 0], MBZ[:, 1], color='b',linewidth=2)
Mkpath_cart = pm.red_to_cart2D(Mrlat3D[0], Mrlat3D[1], Mkpath_red)
ax3.scatter(Mkpath_cart[:,0], Mkpath_cart[:,1])
ax3.arrow(0, 0, Mrlat[0,0], Mrlat[0,1], head_width = 0.01, width = 0.005, color='k', length_includes_head=True)
ax3.arrow(0, 0, Mrlat[1,0], Mrlat[1,1], head_width = 0.01, width = 0.005, color='k', length_includes_head=True)
ax3.axis('off')

plt.show()
fig1.savefig(f"moire_unitcell-{round(angle, 2)}_degree.png", bbox_inches = 'tight', dpi=300)

#---------------------------- plot2: Extracted hopping parameters
fs=12
fig2=plt.figure(figsize=(10, 10))
fig2.subplots_adjust(wspace=0.01,hspace=0.1)
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['font.size'] = fs

ax1 = fig2.add_subplot(321)
ax1.set_aspect('equal', adjustable='box')
daah = ax1.scatter(dA1A2[:,0], dA1A2[:,1], c=tA1A2, cmap='coolwarm', marker="s", s=100)
cbar = plt.colorbar(daah, orientation="vertical",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label="Hopping Parameters (eV)", size=fs)
cbar.ax.tick_params(labelsize=fs)
ax1.set_xlabel(r"dx ($\AA$)", size=fs)
ax1.set_ylabel(r"dy ($\AA$)", size=fs)
ax1.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)


ax2 = fig2.add_subplot(322)
ax2.scatter(dA1A2[:,2],tA1A2,  c='b', marker="s", s=1)
xexA1A2 = np.reshape(np.linspace(0,max(dA1A2[:,0]), 200),(-1,1))
yexA1A2 = np.reshape(np.linspace(0,max(dA1A2[:,1]), 200),(-1,1))
exA1A2 = np.reshape(np.column_stack((xexA1A2, yexA1A2)), (-1,1,2))
rexA1A2 = [np.sqrt(np.dot(x,x)) for x in exA1A2[:,0]]
intp_tA1A2 = pm.matching_interpolate(dA1A2, tA1A2, exA1A2)
ax2.plot(rexA1A2, intp_tA1A2[:,0,0], c='r', lw=2)
ulAA_mag = [np.sqrt(np.dot(x,x)) for x in ulA1A2[:,0,:2]]
ax2.text(((min(ulAA_mag)+max(ulAA_mag))/2)-0.3, 0.16, "Nearest" + "\n" + "neighbor", fontsize = fs)
ax2.axvline(max(ulAA_mag), color='green', lw=1)
ax2.fill_betweenx(np.arange(-0.1, 0.45, 0.01), x1=0, x2=max(ulAA_mag), color='green', alpha=0.2)
uslAA_mag = [np.sqrt(np.dot(x,x)) for x in ulA1A2[:,1,:2]]
ax2.axvline(min(uslAA_mag), color='red', lw=1)
ax2.axvline(max(uslAA_mag), color='red', lw=1)
ax2.fill_betweenx(np.arange(-0.1, 0.45, 0.01), x1=min(uslAA_mag), x2=max(uslAA_mag), color='red', alpha=0.2)
ax2.text(((min(uslAA_mag)+max(uslAA_mag))/2)-0.3, 0.16, "Next nearest" + "\n" + "neighbor", fontsize = fs)
ax2.set_xlim(0.0,3.5)
ax2.set_ylim(-0.1,0.45)
ax2.set_xlabel(r"dr ($\AA$)", size=fs)
ax2.set_ylabel("Hopping Parameters (eV)", size=fs)
ax2.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)

ax3 = fig2.add_subplot(323)
ax3.set_aspect('equal', adjustable='box')
daah = ax3.scatter(dA1B2[:,0], dA1B2[:,1], c=tA1B2, cmap='coolwarm', marker="s", s=100)
cbar = plt.colorbar(daah, orientation="vertical",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label="Hopping Parameters (eV)", size=fs)
cbar.ax.tick_params(labelsize=fs)
ax3.set_xlabel(r"dx ($\AA$)", size=fs)
ax3.set_ylabel(r"dy ($\AA$)", size=fs)
ax3.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)


ax4 = fig2.add_subplot(324)
ax4.scatter(dA1B2[:,2],tA1B2,  c='b', marker="s", s=1)
xexA1B2 = np.reshape(np.linspace(0,max(dA1B2[:,0]), 200),(-1,1))
yexA1B2 = np.reshape(np.linspace(0,max(dA1B2[:,1]), 200),(-1,1))
exA1B2 = np.reshape(np.column_stack((xexA1B2, yexA1B2)), (-1,1,2))
rexA1B2 = [np.sqrt(np.dot(x,x)) for x in exA1B2[:,0]]
intp_tA1B2 = pm.matching_interpolate(dA1B2, tA1B2, exA1B2)
ax4.plot(rexA1B2, intp_tA1B2[:,0,0], c='r', lw=2)
ulAB_mag = [np.sqrt(np.dot(x,x)) for x in ulA1B2[:,0,:2]]
ax4.text(((min(ulAB_mag)+max(ulAB_mag))/2)-0.4, 0.08, 
          "Nearest" + "\n" + "neighbor", fontsize = fs)
ax4.axvline(max(ulAB_mag), color='green', lw=1)
ax4.fill_betweenx(np.arange(-0.1, 0.4, 0.01), 
                  x1=0, x2=max(ulAB_mag), color='green', alpha=0.2)
uslAB_mag = [np.sqrt(np.dot(x,x)) for x in ulA1B2[:,1,:2]]
ax4.text(((min(uslAB_mag)+max(uslAB_mag))/2)-0.3, 0.20, 
          "Next nearest" + "\n" + "neighbor", fontsize = fs)
ax4.axvline(min(uslAB_mag), color='red', lw=1)
ax4.axvline(max(uslAB_mag), color='red', lw=1)
ax4.fill_betweenx(np.arange(-0.1, 0.4, 0.01), 
                  x1=min(uslAB_mag), x2=max(uslAB_mag), color='red', alpha=0.2)
ax4.set_xlim(0.0,3.5)
ax4.set_ylim(-0.1,0.4)
ax4.set_xlabel(r"dr ($\AA$)", size=fs)
ax4.set_ylabel("Hopping Parameters (eV)", size=fs)
ax4.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)

ax5 = fig2.add_subplot(325)
ax5.set_aspect('equal', adjustable='box')
daah = ax5.scatter(dB1A2[:,0], dB1A2[:,1], c=tB1A2, cmap='coolwarm', marker="s", s=100)
cbar = plt.colorbar(daah, orientation="vertical",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label="Hopping Parameters (eV)", size=fs)
cbar.ax.tick_params(labelsize=fs)
ax5.set_xlabel(r"dx ($\AA$)", size=fs)
ax5.set_ylabel(r"dy ($\AA$)", size=fs)
ax5.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)


ax6 = fig2.add_subplot(326)
ax6.scatter(dB1A2[:,2],tB1A2,  c='b', marker="s", s=1)
xexB1A2 = np.reshape(np.linspace(0,max(dB1A2[:,0]), 200),(-1,1))
yexB1A2 = np.reshape(np.linspace(0,max(dB1A2[:,1]), 200),(-1,1))
exB1A2 = np.reshape(np.column_stack((xexB1A2, yexB1A2)), (-1,1,2))
rexB1A2 = [np.sqrt(np.dot(x,x)) for x in exB1A2[:,0]]
intp_tB1A2 = pm.matching_interpolate(dB1A2, tB1A2, exB1A2)
ax6.plot(rexB1A2, intp_tB1A2[:,0,0], c='r', lw=2)
ulBA_mag = [np.sqrt(np.dot(x,x)) for x in ulB1A2[:,0,:2]]
ax6.text(((min(ulBA_mag)+max(ulBA_mag))/2)-0.4, 0.08, 
          "Nearest" + "\n" + "neighbor", fontsize = fs)
ax6.axvline(max(ulBA_mag), color='green', lw=1)
ax6.fill_betweenx(np.arange(-0.1, 0.4, 0.01), 
                  x1=0, x2=max(ulBA_mag), color='green', alpha=0.2)
uslBA_mag = [np.sqrt(np.dot(x,x)) for x in ulB1A2[:,1,:2]]
ax6.text(((min(uslBA_mag)+max(uslBA_mag))/2)-0.3, 0.20, 
          "Next nearest" + "\n" + "neighbor", fontsize = fs)
ax6.axvline(min(uslBA_mag), color='red', lw=1)
ax6.axvline(max(uslBA_mag), color='red', lw=1)
ax6.fill_betweenx(np.arange(-0.1, 0.4, 0.01), 
                  x1=min(uslBA_mag), x2=max(uslBA_mag), color='red', alpha=0.2)
ax6.set_xlim(0.0,3.5)
ax6.set_ylim(-0.1,0.4)
ax6.set_xlabel(r"dr ($\AA$)", size=fs)
ax6.set_ylabel("Hopping Parameters (eV)", size=fs)
ax6.tick_params(axis ='both', which ='major', length = 4, direction='in', labelsize=fs)


fig2.tight_layout()
plt.show()
fig2.savefig(f"Hopping_analysis_{round(angle, 2)}_degree.png", dpi=300)

# #---------------------------- plot3: Total-energy, VdW-energy, Fermi-energy and Bandgap-energy on moire structure
fs = 60
fig3 = plt.figure(figsize=(15, 10))
fig3.subplots_adjust(wspace=0.1, hspace=0.0)
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = fs

ax1 = fig3.add_subplot(141)
ax1.set_aspect('equal', adjustable='box')
Emap= plt.scatter(RA[:,0], RA[:,1], c=MTE[:,0], cmap='coolwarm', marker="s", s=fs, alpha=1)

cbar = plt.colorbar(Emap, orientation="horizontal",extend = 'neither',
                    pad=0.01, shrink=0.9, aspect=20, format='%.1f', 
                    ticks=[((max(MTE[:,0])+min(MTE[:,0]))/2+min(MTE[:,0]))/2, ((max(MTE[:,0])+min(MTE[:,0]))/2+max(MTE[:,0]))/2])
cbar.set_label(label="Local Total Energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax1.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
ax1.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax1.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax1.set_xticks([])
ax1.set_yticks([])

ax2 = fig3.add_subplot(142)
ax2.set_aspect('equal', adjustable='box')
Emap= plt.scatter(RA[:,0], RA[:,1], c=MVdW[:,0], cmap='coolwarm', marker="s", s=fs, alpha=1)
cbar = plt.colorbar(Emap, orientation="horizontal",extend = 'neither',
                    pad=0.01, shrink=0.9, aspect=20, format='%.0f')
cbar.set_label(label="Local VdW Energy (meV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax2.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
ax2.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax2.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax2.set_xticks([])
ax2.set_yticks([])


ax3 = fig3.add_subplot(143)
ax3.set_aspect('equal', adjustable='box')
Emap= plt.scatter(RA[:,0], RA[:,1], c=MFE[:,0], cmap='coolwarm', marker="s", s=fs, alpha=1)
cbar = plt.colorbar(Emap, orientation="horizontal",extend = 'neither',
                    pad=0.01, shrink=0.9, aspect=20, format='%.2f')
cbar.set_label(label="Local Fermi Energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax3.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
ax3.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax3.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax3.set_xticks([])
ax3.set_yticks([])

ax4 = fig3.add_subplot(144)
ax4.set_aspect('equal', adjustable='box')
Emap= plt.scatter(RA[:,0], RA[:,1], c=Mbg[:,0], cmap='coolwarm', marker="s", s=fs, alpha=1)
cbar = plt.colorbar(Emap, orientation="horizontal",extend = 'neither',
                    pad=0.01, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label="Local Band Gap Energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax4.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
ax4.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax4.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax4.set_xticks([])
ax4.set_yticks([])

plt.show()
fig3.savefig(f"total-VdW-Fermi_Energy_{round(angle, 2)}_degree.png", bbox_inches = 'tight', dpi=300)

#----------------------------plot4: Mapped hopping parameters on moire structure for each orbital coupling in two neighbors
fs=60
fig4=plt.figure(figsize=(18, 12))
fig4.subplots_adjust(wspace=0.0,hspace=0.0)
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = fs

ax1 = fig4.add_subplot(3,4,1)
ax1.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RA[:,0], RA[:,1], c=MtA1A2[:,0], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label=r"nn AA$^\prime$ Coupling Energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax1.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RA)):
    ax1.arrow(RA[i,0], RA[i,1],lA1A2[i,0,0],lA1A2[i,0,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax1.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax1.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax1.set_xticks([])
ax1.set_yticks([])

ax2 = fig4.add_subplot(3,4,2)
ax2.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RA[:,0], RA[:,1], c=MtA1B2[:,0], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label=r"nn AB$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax2.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RA)):
    ax2.arrow(RA[i,0], RA[i,1],lA1B2[i,0,0],lA1B2[i,0,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax2.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax2.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax2.set_xticks([])
ax2.set_yticks([])

ax3 = fig4.add_subplot(3,4,3)
ax3.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RB[:,0], RB[:,1], c=MtB1B2[:,0], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label=r"nn BB$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax3.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RB)):
    ax3.arrow(RB[i,0], RB[i,1],lB1B2[i,0,0],lB1B2[i,0,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax3.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax3.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax3.set_xticks([])
ax3.set_yticks([])

ax4 = fig4.add_subplot(3,4,4)
ax4.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RB[:,0], RB[:,1], c=MtB1A2[:,0], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.1f')
cbar.set_label(label=r"nn BA$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax4.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RB)):
    ax4.arrow(RB[i,0], RB[i,1],lB1A2[i,0,0],lB1A2[i,0,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax4.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax4.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax4.set_xticks([])
ax4.set_yticks([])

ax5 = fig4.add_subplot(3,4,5)
ax5.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RA[:,0], RA[:,1], c=MtA1A2[:,1], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.2f')
cbar.set_label(label=r"nnn AA$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax5.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RA)):
    ax5.arrow(RA[i,0], RA[i,1],lA1A2[i,1,0],lA1A2[i,1,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax5.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax5.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax5.set_xticks([])
ax5.set_yticks([])

ax6 = fig4.add_subplot(3,4,6)
ax6.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RA[:,0],RA[:,1], c=MtA1B2[:,1], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.2f')
cbar.set_label(label=r"nnn AB$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax6.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RA)):
    ax6.arrow(RA[i,0], RA[i,1],lA1B2[i,1,0],lA1B2[i,1,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax6.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax6.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax6.set_xticks([])
ax6.set_yticks([])

ax7 = fig4.add_subplot(3,4,7)
ax7.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RB[:,0], RB[:,1], c=MtB1B2[:,1], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.2f')
cbar.set_label(label=r"nnn BB$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax7.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RB)):
    ax7.arrow(RB[i,0], RB[i,1],lB1B2[i,1,0],lB1B2[i,1,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax7.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax7.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax7.set_xticks([])
ax7.set_yticks([])

ax8 = fig4.add_subplot(3,4,8)
ax8.set_aspect('equal', adjustable='box')
allatm = plt.scatter(RB[:,0], RB[:,1], c=MtB1A2[:,1], cmap='coolwarm', marker="s", s=fs)
cbar = plt.colorbar(allatm, orientation="horizontal",extend = 'neither',
                    pad=0.02, shrink=0.9, aspect=20, format='%.2f')
cbar.set_label(label=r"nnn BA$^\prime$ Coupling energy (eV)", size=15)
cbar.ax.tick_params(labelsize=15)
ax8.plot(Muc[:,0],Muc[:,1], color='k',linewidth=1)
for i in range(len(RB)):
    ax8.arrow(RB[i,0], RB[i,1],lB1A2[i,1,0],lB1A2[i,1,1],head_width = 0.3, width = 0.06, color='k', alpha= 1, length_includes_head=True)
ax8.set_xlim(xmin=-Muc[2,0]*0.1,xmax=Muc[2,0]*1.1)
ax8.set_ylim(ymin=-Muc[2,1]*0.1,ymax=Muc[2,1]*1.1)
ax8.set_xticks([])
ax8.set_yticks([])

fig4.tight_layout()
plt.show()
fig4.savefig(f"nearest_and_next_nearest_hoppings_{round(angle, 2)}_degree.png", bbox_inches = 'tight', dpi=300)

#----------------------------plot5: Band-structure and DOS plot
k = list(range(E_Ef.shape[1]))

fs=10
fig5=plt.figure(figsize=(14, 8))
fig5.subplots_adjust(wspace=0.0,hspace=0.0)
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['font.size'] = fs
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.serif"] = "Times New Roman"
ax1 = fig5.add_subplot(121)
ax1.set_xlim(xmin=0,xmax=hs[-1])
ax1.set_ylim(ymin=-1,ymax=1)
ax1.set_ylabel('$E-E_{Fermi}$ (eV)',fontsize=fs+20)
ax1.set_xlabel('kpoints',fontsize=fs+20)
ax1.axhline(y=0.0, color='b', linestyle='--')
for vl in hs:
  ax1.axvline(vl, color='black', lw=1)
for i in range(E_Ef.shape[0]):
  ax1.plot(k, E_Ef[i],"k-",zorder=-50, alpha=1)

ax1.set_xticks(hs)
ax1.set_xticklabels(hstitle)
ax1.tick_params(axis ='both', which ='major', 
                length = 4, direction='in', labelsize=30)
ax1.plot(list(range(E_Ef.shape[1])), E_Ef[int(dim/2)-2,:],"r-")

# ax2 = fig5.add_subplot(122)
# ax2.plot(dos, energy-Ef[round(angle,2)], "k-")
# ax2.axhline(y=0.0, color='b', linestyle='--')
# ax2.set_ylim(-0.2, 0.2)
# ax2.set_xlim(0, max(dos)+100)
# ax2.set_xlabel("DOS",fontsize=fs+20)
# plt.xticks([])
# plt.yticks([])

fig5.tight_layout()
plt.show()
fig5.savefig(f"Moire_band_{round(angle, 2)}-degree.png", bbox_inches = 'tight', dpi=300)
