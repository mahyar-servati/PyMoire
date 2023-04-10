#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:25:37 2022

@author: Mahyar Servati
"""
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate.rbf import Rbf

pi = np.pi

#--------------------------------------------Unit-cell
def unit_cell(lat):
    Ruc = np.array([np.array([0.0, 0.0]), lat[0], lat[0]+lat[1], lat[1], np.array([0.0, 0.0])])
    return Ruc
#--------------------------------------------Real lattice to reciprocal lattice
def real_to_reciprocal(lat):
    rec_lat = np.zeros_like(lat)
    V = np.linalg.det(lat)
    rec_lat[0] = 2*pi*np.array([lat[1][1],(-1)*lat[1][0]])/V
    rec_lat[1] = 2*pi*np.array([(-1)*lat[0][1],lat[0][0]])/V
    return rec_lat
#--------------------------------------------First Briloen Zone
def first_brillouin(rlat):
    BZ = np.zeros((6+1,2))
    ra = np.sqrt(np.dot(rlat[0],rlat[0]))
    rb = np.sqrt(np.dot(rlat[1],rlat[1]))
    alpha = np.degrees(np.arccos(np.dot(rlat[0],rlat[1])/(ra*rb)))
    if alpha < 90:
        K1x = (rlat[0][0] + rlat[1][0])/3
        K1y = (rlat[0][1] + rlat[1][1])/3
        for i in range(6+1):
           BZ[i][0] = K1x*np.cos(np.radians(60*i))-K1y*np.sin(np.radians(60*i))
           BZ[i][1] = K1x*np.sin(np.radians(60*i))+K1y*np.cos(np.radians(60*i))
    else:
        K1x = (2*rlat[0][0] + rlat[1][0])/3
        K1y = (2*rlat[0][1] + rlat[1][1])/3
        for i in range(6+1):
           BZ[i][0] = K1x*np.cos(np.radians(60*i))-K1y*np.sin(np.radians(60*i))
           BZ[i][1] = K1x*np.sin(np.radians(60*i))+K1y*np.cos(np.radians(60*i))
    return BZ
#----------------------------------------------------Reference Atom positions
def ref_latt(lat, atm, begin, end):
    RA=np.zeros(((begin+end)*(begin+end),2))
    RB=np.zeros(((begin+end)*(begin+end),2))
    sa = atm[0]
    sb = atm[1]
    n = 0
    for i in np.arange(-1*begin, end):
        for j in np.arange(-1*begin, end):
            RA[n] = np.array([j*lat[0][0]+i*lat[1][0]+sa[0], j*lat[0][1]+i*lat[1][1]+sa[1]])
            RB[n] = np.array([j*lat[0][0]+i*lat[1][0]+sb[0], j*lat[0][1]+i*lat[1][1]+sb[1]])
            n +=1
    return RA, RB

def trunc_latt(MA, RA, RB):
    trunc_RA=[]
    trunc_RB=[]
    for i in range(len(RA)):
        if (RA[i,1] > -2*2.46) and (RA[i,1] < MA+0.1*MA):
            trunc_RA.append(RA[i])
            trunc_RB.append(RB[i])
    trunc_RA = np.array(trunc_RA)
    trunc_RB = np.array(trunc_RB)
    return trunc_RA, trunc_RB

#----------------------------------------------------Twist layers

def twist(teta, lat):
    tlat = np.zeros_like(lat)
    rotation = np.array([[np.cos(teta), np.sin(teta)],
                        [-np.sin(teta), np.cos(teta)]])
    for i in range(len(lat)):
        tlat[i] = np.dot(lat[i], rotation)
    return tlat

def twist3D(teta, lat):
    tlat = np.zeros_like(lat)
    rotation = np.array([[np.cos(teta), np.sin(teta), 0.0],
                        [-np.sin(teta), np.cos(teta), 0.0],
                        [0.00000000000, 0.0000000000, 1.0]])
    for i in range(len(lat)):
        tlat[i] = np.dot(lat[i], rotation)
    return tlat
#----------------------------------------------------Lattice
def Lattice(n):
    L = []
    for i in np.arange(-n, n+1):
        for j in np.arange(-n, n+1):
            L.append([i, j])
    return L
#----------------------------------------------------Moire Lattice
#--------------------------------------------Commensurate Rotation
def comensorate_moire(x):
    angle = np.degrees(np.arccos((3*(x*x) + 3*x + 0.5)/(3*(x*x) + 3*x + 1)))
    return angle
#--------------------------------------------Incommensurate Rotation
def incomensorate_moire(angle, A):
    teta = np.radians(angle)
    M_angle = angle/2
    MA = (1/(2*np.sin(teta/2)))*A
    return teta, MA, M_angle
#--------------------------------------------Moire vectors
def moire_vec(MA, M_angle):
    Mlat = np.zeros((2,2))
    Mlat[0] = np.array([MA*np.cos(np.radians(M_angle+60)),
                       MA*np.sin(np.radians(M_angle+60))])
    Mlat[1] = np.array([MA*np.cos(np.radians(M_angle)),
                        MA*np.sin(np.radians(M_angle))])
    return Mlat

 #--------------------------------------------k-path
def k_path(kinfo):
    kinfo = np.array(kinfo)
    kvec = kinfo[:,:3].astype(float)
    kpath = [kvec[0]]
    points = kinfo[:,4].astype(int)
    hs = [0]
    p = 0
    for k in range(1,len(kvec)):
        klist = np.linspace(kvec[k-1], kvec[k], points[k-1])
        for l in range(len(klist)):
            if (np.array(kpath[-1]) != klist[l]).any():
                kpath.append(klist[l])
        p += points[k-1]
        hs.append(p-1)
    hstitle = kinfo[:,3]
    return np.array(kpath), hs, hstitle


#---------------------------------------------------- Interlayer distance between non-twisted and twisted atoms

def distance(A, RA, TB): 
    neighbor_num = 3
    lab=np.zeros((len(RA),neighbor_num, 3), dtype=float)
    for i in range(len(RA)):
        dis=[]
        RAr = np.sqrt(np.dot(RA[i],RA[i]))
        for j in range(len(TB)):
            TBr = np.sqrt(np.dot(TB[j],TB[j]))
            if (abs(TBr-RAr) < 2*A):
                disx = TB[j][0]-RA[i][0]
                disy = TB[j][1]-RA[i][1]
                disr = np.sqrt(pow(disx,2)+pow(disy,2))
                dis.append([disx, disy, disr])
        sort_ab = sorted(dis, key=lambda x: x[2])
        for n in range(neighbor_num):
            lab[i,n] = np.array(sort_ab[n])   
    return lab

def intra_distance(Rlat, Ratm):
    Lat = Lattice(1)
    if (Ratm[0] == 0.0 and Ratm[1] == 0.0):
        rA = np.sqrt(np.dot(Rlat[0],Rlat[0]))
    else:
        rA = np.sqrt(np.dot(Ratm,Ratm))
    nighbors = []

    for l in Lat:
        vec = Ratm+(l[0]*Rlat[0])+(l[1]*Rlat[1])
        if 0.0 < np.sqrt(np.dot(vec,vec)) < rA+0.01:
            nighbors.append(vec)
    return np.array(nighbors)

#---------------------------------------------------Moire-ribbon

def ribbon(RA, la, Mlat, l):
    rb_la = []
    red = cart_to_red_2D(Mlat, RA)  
    
    for i in range(len(RA)):
        if red[i,1]-l <= red[i,0] <= red[i,1]+l and red[i,0]-l <= red[i,1] <= red[i,0]+l:
            rb_la.append(la[i])
    rb_la = np.array(rb_la)
    return rb_la

#----------------------------------------------------Moire Unitcell

def moire_unitcell(RA, la, Mlat):
    ula=[]
    red = cart_to_red_2D(Mlat, RA)    
    for i in range(len(RA)):
        if 0 <= red[i,0] < 1 and 0 <= red[i,1] < 1:
            ula.append(la[i])
    ula = np.array(ula)
    return ula

#--------------------------------------------Extracting hopping parameters
def extract_hopping(B, i, file_direction, column):
    with open(file_direction, "r") as f:
        lines= f.readlines()
    dx=np.zeros((len(lines),1))
    dy=np.zeros((len(lines),1))
    dr=[]
    t=np.zeros((len(lines),1))
    for n in range(len(lines)):
        dx[n]=(B+float(lines[n].split()[0]))
        dy[n]=float(lines[n].split()[1])
        t[n]=(i)*float(lines[n].split()[column])
    
    for i in range(len(dx)):
        dr.append(np.sqrt(pow(dx[i], 2) + pow(dy[i], 2)))
    d = np.column_stack((dx, dy, dr))
    return d, t

def intra_hopping(file_direction):
    with open(file_direction, "r") as f:
        lines= f.readlines()
    t = [round(float(i),4) for i in lines[0].split()]
    t.pop(0)
    t.pop(0)
    tset = np.array(list(set(t)))
    tabs = np.array([abs(i) for i in tset])
    targsort = np.argsort(tabs)
    tsort = tset[targsort][::-1]
    return tsort

#----------------------------------------------------Outliers fixing
def outlier(dB1A2, tB1A2):
    knn = 1.2
    nnrange = 0.08

    for i in range(len(dB1A2)):
        nnt = []
        for j in range(len(dB1A2)):
            if (abs(dB1A2[i,2]-dB1A2[j,2]) > 0.0 and abs(dB1A2[i,2]-dB1A2[j,2]) < knn*abs(dB1A2[0,2]-dB1A2[1,2])):
                nnt.append(tB1A2[j,0])
        nnmean = np.mean(nnt)
        # nnstd = np.std(nnt)
        if (tB1A2[i,0] > nnmean+nnrange or tB1A2[i,0] < nnmean-nnrange):
            tB1A2[i,0] = nnmean
    return tB1A2
#----------------------------------------------------Matching Hopping parameters into moire lattice
def matching_interpolate(dA1A2, tA1A2, lAB):
    rA1A2_1d = []; tA1A2_1d =[]
    tAB = np.zeros((len(lAB),len(lAB[0]), 1), dtype=float)
    for i in range(len(dA1A2)):
        if (dA1A2[i,0] == dA1A2[i,1] and dA1A2[i,0]>=0 and dA1A2[i,1]>=0):
            rA1A2_1d.append(dA1A2[i,2])
            tA1A2_1d.append(tA1A2[i,0])
    y_f = interp1d(rA1A2_1d, tA1A2_1d, 'cubic') 
    rlAB = np.zeros((len(lAB),len(lAB[0]), 1), dtype=float)
    for i in range(len(lAB)):
        for l in range(len(lAB[i])):
            rlAB[i,l] = np.sqrt(pow(lAB[i,l,0],2)+pow(lAB[i,l,1],2))
            if rlAB[i,l] > max(rA1A2_1d):
                rlAB[i,l] = max(rA1A2_1d)
            tAB[i,l] = y_f(rlAB[i,l])
    return tAB

def matching_interpolate2D(dA1B2, tA1B2, lAB):
    dx = list(set(dA1B2[:,0]))
    dx.sort()
    dy = dx
    h = np.reshape(tA1B2, (-1, len(dx))).T
    y_f = interp2d(dx, dy, h, kind='cubic') 
    t = np.zeros((len(lAB),len(lAB[0]),1), dtype=float)
    for i in range(len(lAB)):
        for l in range(len(lAB[i])):
            rlAB = np.sqrt(pow(lAB[i,l,0],2)+pow(lAB[i,l,1],2))
            if rlAB < max(dA1B2[:,2]):
                t[i,l] = y_f(abs(lAB[i,l,0]), abs(lAB[i,l,1]))
            else:
                t[i,l] = np.array([np.mean(tA1B2)])
    return t

def energy_interpolate2D(dTE, TE, lAB):
    dx = list(set(dTE[:,0]))
    dx.sort()
    dy = dx
    h = np.reshape(TE, (-1, len(dx))).T
    y_f = interp2d(dx, dy, h, kind='cubic') 
    t = np.zeros((len(lAB),1), dtype=float)
    for i in range(len(lAB)):
        rlAB = np.sqrt(pow(lAB[i,0,0],2)+pow(lAB[i,0,1],2))
        if rlAB < max(dTE[:,2]):
            t[i] = y_f(abs(lAB[i,0,0]), abs(lAB[i,0,1]))
        else:
            t[i] = np.array([np.mean(TE)])
    return t

def matching_Rbf(dA1B2, tA1B2, lAB):
    x = dA1B2[:,0]
    y = dA1B2[:,1]
    y_f = Rbf(x, y, tA1B2, function='multiquadric')
    t =np.zeros((len(lAB), len(lAB[0]), 1), dtype=float)
    for i in range(len(lAB)):
        for l in range(len(lAB[i])):
            rlAB = np.sqrt(pow(lAB[i,l,0],2)+pow(lAB[i,l,1],2))
            if rlAB < max(dA1B2[:,2]):
                t[i,l] = y_f(lAB[i,l,0], lAB[i,l,1])
            else:
                t[i,l] = np.array([np.mean(tA1B2)])
    return t
#--------------------------------------------Change Coordinates

def cart_to_red(lat, cart):
    cnv=np.array(lat)
    cnv=cnv.T
    cnv=np.linalg.inv(cnv)
    red=np.zeros_like(cart,dtype=float)
    for i in range(0,len(cart)):
        for l in range(len(cart[i])):
            red[i,l]=np.dot(cnv,cart[i,l])
    return red

def cart_to_red_2D(lat, cart):
    cnv=np.array(lat)
    cnv=cnv.T
    cnv=np.linalg.inv(cnv)
    red=np.zeros_like(cart,dtype=float)
    for i in range(0,len(cart)):
        red[i]=np.dot(cnv,cart[i])
    return red

def red_to_cart(tmp,red):
    (a1,a2,a3)=tmp
    cart=np.zeros_like(red,dtype=float)
    for i in range(0,len(cart)):
        cart[i,:]=a1*red[i][0]+a2*red[i][1]+a3*red[i][2]
    return cart

def red_to_cart2D(a1, a2, red):
    cart=np.zeros_like(red, dtype=float)
    for i in range(0,len(cart)):
        cart[i,:]=a1*red[i][0]+a2*red[i][1]
    return cart

#---------------------------------------------------- Onsite
def onsite(red_uRA, eA1A1, norb, kind='None'):
    if  kind.upper() == "A1":
        index1 = 0
    if  kind.upper() == "B1":
        index1 = norb[0]
    if  kind.upper() == "A2":
        index1 = norb[0]+norb[1]
    if  kind.upper() == "B2":
        index1 = norb[0]+norb[1]+norb[2]
    onsite_index = []
    for i, R in enumerate(red_uRA):
        onsite_index.append([0, 0, 0, index1+i, index1+i, eA1A1, np.array([0.0, 0.0, 0.0])])
    return onsite_index    

    
#---------------------------------------------------- coupled orbitals indexing

def hopping_index(red_uRA, red_uRB, red_distances, norb, kind='None'):
    if  kind[:2].upper() == "A1":
        index1 = 0
    if  kind[:2].upper() == "B1":
        index1 = norb[0]
    if  kind[:2].upper() == "A2":
        index1 = norb[0]+norb[1]
    if  kind[:2].upper() == "B2":
        index1 = norb[0]+norb[1]+norb[2]
        
    if  kind[-2:].upper() == "A1":
        index2 = 0
    if  kind[-2:].upper() == "B1":
        index2 = norb[0]
    if  kind[-2:].upper() == "A2":
        index2 = norb[0]+norb[1]
    if  kind[-2:].upper() == "B2":
        index2 = norb[0]+norb[1]+norb[2]
      
    if kind[1]==kind[-1]:
        ws = np.zeros((len(red_uRA),3*len(red_distances)), dtype=int)
        hop_index =np.zeros((len(red_uRA),1+len(red_distances)), dtype=int)
    if not kind[1]==kind[-1]:
        ws = np.zeros((len(red_uRA),3), dtype=int)
        hop_index =np.zeros((len(red_uRA),2), dtype=int)
    for i, R in enumerate(red_uRA):
        hop_index[i,0]=i+index1
        if kind[1]==kind[-1]:
            red_inter_atom_dis = red_distances
        if not kind[1]==kind[-1]:
            red_inter_atom_dis = np.array([red_distances[i]])
        for j, t in enumerate(red_inter_atom_dis):
            redd = R+t
            if redd[0]>1.0:
                redd[0] = redd[0]-1
                ws[i,j*3]=1
            if redd[0]<0.0:
                redd[0] = redd[0]+1
                ws[i,j*3]=-1
            if redd[1]>1.0:
                redd[1] = redd[1]-1
                ws[i,j*3+1]=1
            if redd[1]<0.0:
                redd[1] = redd[1]+1
                ws[i,j*3+1]=-1
            for k, Rb in enumerate(red_uRB):
                if (abs(Rb[0]-redd[0])<1.0E-8) and (abs(Rb[1]-redd[1])<1.0E-8):
                    hop_index[i,j+1]= k+index2
            # prevent index from remaining zero
            if sum(norb) > 10000:
                if hop_index[i,j+1]<index2 or hop_index[i,j+1]>index2+len(red_uRB):
                    hop_index[i,j+1] = index2
    return hop_index, ws
#----------------------------------------------------Stacking data
def stack(ws, hi, hop, red_len):
    hoppings = []
    if len(red_len) == len(hi[0])-1:
        hop = hop*np.ones((len(hi)))
        for i in range(len(hi[0])-1):
            for j in range(len(hi)):
                hoppings.append([ws[j,3*i], ws[j,3*i+1], ws[j,3*i+2], hi[j,0], hi[j,i+1], hop[j], red_len[i]])
    if not len(red_len) == len(hi[0])-1:
        for i in range(len(hi)):
            hoppings.append([ws[i,0], ws[i,1], ws[i,2], hi[i,0], hi[i,1], hop[i,0], red_len[i]])
            # hoppings.append([0, 0, 0, hi[i,1], hi[i,0], hop[i,0], -1*red_len[i]])
            
    return hoppings
#----------------------------------------------------Save to file
def save_to_file(direction, lists):
    list1d = lists.reshape(-1)
    file = open(direction, "w")
    for li in list1d:
        file.write(str(li) + " ")
    file.close()
#----------------------------------------------------Open from file
def open_from_file(direction, lists):
    with open(direction, "r") as file:
        lines= file.readlines()    
    list1d = np.array(list(map(float, lines[0].strip().split())))
    flist = np.reshape(list1d, (len(lists), 3, 3))
    return flist

#----------------------------------------------------hr file
def neigh_extract(red_inter_AA):
    des = np.zeros((len(red_inter_AA), len(red_inter_AA[0,0])), dtype=float)
    sdes = np.zeros((len(red_inter_AA), len(red_inter_AA[0,0])), dtype=float)
    ssdes = np.zeros((len(red_inter_AA), len(red_inter_AA[0,0])), dtype=float)
    for i in range(len(red_inter_AA)):
        des[i] = red_inter_AA[i,0]
        sdes[i] = red_inter_AA[i,1]
        ssdes[i] = red_inter_AA[i,2]
    return des, sdes, ssdes
#----------------------------------------------------
def progress_bar(current, total, bar_length=20):
    fraction = (current+1) / total

    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '

    ending = '\n' if current+1 == total else ''

    print('\r', f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)
