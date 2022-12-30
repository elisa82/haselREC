import numpy as np
import pandas as pd
import os
from .scale_acc import scale_acc
from .create_output_files import create_output_files
from zipfile import ZipFile
from os.path import basename

class NTC08_param:
  def __init__(self, lon, lat, ag, f0, Tcstar):
    self.lon = lon
    self.lat = lat
    self.ag = ag
    self.f0 = f0
    self.Tcstar = Tcstar

def select_NTC08(sitelon,sitelat,ag,f0,Tcstar,longnodo,latnodo):
    dist=[]
    for i in range(len(longnodo)):
        #Viene preso il punto delle NTC piu' vicino
        dist.append(np.sqrt((sitelat-latnodo[i])**2+(sitelon-longnodo[i])**2))
    dist=np.array(dist)
    index_min = np.argmin(dist)
    ag_site = ag[index_min]
    f0_site = f0[index_min]
    Tcstar_site = Tcstar[index_min]
    return ag_site,f0_site,Tcstar_site

def compute_NTC08_design_spectrum(ag,f0,Tcstar,x):
    ag=ag/10.
    Smorzamento=5.
    Cc=1
    ST=1 #amplificazione topografica
    Ss=1 #amplificazione stratigrafica
    S=Ss*ST #e' il coefficiente che tiene conto della categoria di sottosuolo e delle condizioni topografiche
    TC=Cc*Tcstar #e' il periodo corrispondente al'inizio del tratto a velocita' costante dello spettro
    TD=ag*4.+1.6 #e' il periodo corrispondente al'inizio del tratto a spostamento costante dello spettro
    TB=TC/3. #e' il periodo corrispondente al'inizio del tratto dello spettro ad accelerazione costante
    #eta e' il fattore che altera lo spettro elastico per coefficienti di smorzamento viscosi convenzionali diversi dal 5%
    EtaTENT=np.sqrt(10./(5.+Smorzamento))
    if EtaTENT>0.55:
        Eta=EtaTENT
    else:
        Eta=0.55
    npoints=len(x)
    y=np.zeros((npoints))
    for i in range(npoints):
        if(x[i]<TB):
            y[i]=ag*S*Eta*f0*(x[i]/TB+(1.-x[i]/TB)/(Eta*f0))
        if(x[i]>=TB and x[i]<TC):
            y[i]=ag*S*Eta*f0
        if x[i]>=TC and x[i]<TD:
            y[i]=ag*S*Eta*f0*(TC/x[i])
        if x[i]>=TD:
            y[i]=ag*S*Eta*f0*(TC*TD/(x[i])**2)
    return y


def read_allegatoB(file_input):
    allegatoB=pd.read_csv(file_input,header=None,skiprows=1,names=('lon','lat', 'ag_TR30', 'f0_TR30', 'Tcstar_TR30', 'ag_TR50', 'f0_TR50', 'Tcstar_TR50', 'ag_TR72', 'f0_TR72', 'Tcstar_TR72', 'ag_TR101', 'f0_TR101', 'Tcstar_TR101', 'ag_TR140', 'f0_TR140', 'Tcstar_TR140', 'ag_TR201', 'f0_TR201', 'Tcstar_TR201', 'ag_TR475', 'f0_TR475', 'Tcstar_TR475', 'ag_TR975', 'f0_TR975', 'Tcstar_TR975', 'ag_TR2475', 'f0_TR2475', 'Tcstar_TR2475'),delimiter='\t')	
    ag = np.zeros((len(allegatoB),9))
    f0 = np.zeros((len(allegatoB),9))
    Tcstar = np.zeros((len(allegatoB),9)) 
    latnodo = np.zeros((len(allegatoB))) 
    longnodo = np.zeros((len(allegatoB))) 
    for i in range(len(allegatoB)):
        ag[i,0] = allegatoB['ag_TR30'][i]
        ag[i,1] = allegatoB['ag_TR50'][i]
        ag[i,2] = allegatoB['ag_TR72'][i]
        ag[i,3] = allegatoB['ag_TR101'][i]
        ag[i,4] = allegatoB['ag_TR140'][i]
        ag[i,5] = allegatoB['ag_TR201'][i]
        ag[i,6] = allegatoB['ag_TR475'][i]
        ag[i,7] = allegatoB['ag_TR975'][i]
        ag[i,8] = allegatoB['ag_TR2475'][i]
        f0[i,0] = allegatoB['f0_TR30'][i]
        f0[i,1] = allegatoB['f0_TR50'][i]
        f0[i,2] = allegatoB['f0_TR72'][i]
        f0[i,3] = allegatoB['f0_TR101'][i]
        f0[i,4] = allegatoB['f0_TR140'][i]
        f0[i,5] = allegatoB['f0_TR201'][i]
        f0[i,6] = allegatoB['f0_TR475'][i]
        f0[i,7] = allegatoB['f0_TR975'][i]
        f0[i,8] = allegatoB['f0_TR2475'][i]
        Tcstar[i,0] = allegatoB['Tcstar_TR30'][i]
        Tcstar[i,1] = allegatoB['Tcstar_TR50'][i]
        Tcstar[i,2] = allegatoB['Tcstar_TR72'][i]
        Tcstar[i,3] = allegatoB['Tcstar_TR101'][i]
        Tcstar[i,4] = allegatoB['Tcstar_TR140'][i]
        Tcstar[i,5] = allegatoB['Tcstar_TR201'][i]
        Tcstar[i,6] = allegatoB['Tcstar_TR475'][i]
        Tcstar[i,7] = allegatoB['Tcstar_TR975'][i]
        Tcstar[i,8] = allegatoB['Tcstar_TR2475'][i]
        latnodo[i] = allegatoB['lat'][i]
        longnodo[i] = allegatoB['lon'][i]
    node = NTC08_param(longnodo, latnodo, ag, f0, Tcstar)
    return(node)

def create_NTC18(sitelon,sitelat,TR,x,file_AllegatoB):
    AllegatoB=read_allegatoB(file_AllegatoB)
    if TR==30:
        pos=0
    if TR==50:
        pos=1
    if TR==72:
        pos=2
    if TR==101:
        pos=3
    if TR==140:
        pos=4
    if TR==201:
        pos=5
    if TR==475:
        pos=6
    if TR==975:
        pos=7
    if TR==2475:
        pos=8
    ag,f0,Tcstar = select_NTC08(sitelon,sitelat,AllegatoB.ag[:,pos],AllegatoB.f0[:,pos],AllegatoB.Tcstar[:,pos],AllegatoB.lon,AllegatoB.lat)
    NTC08_spectrum = compute_NTC08_design_spectrum(ag,f0,Tcstar,x)
    return NTC08_spectrum
