import os
import _thread
from scipy.interpolate import interp1d
from bs4 import BeautifulSoup
import requests
from openpyxl import load_workbook
from openpyxl import Workbook
from mendeleev import element
from numpy import loadtxt, log10, array
import numpy as np
import scipy
from scipy import constants,optimize
from scipy.interpolate import InterpolatedUnivariateSpline
from chempy import Substance
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd
from pandas import read_csv
import math
import time
import sqlalchemy
import eel
import json
from datetime import datetime
import sys
sys.setrecursionlimit(5000)
n = scipy.constants.N_A

class Compound:
    def __init__(self, c_en, ic_mat, den_mat, mfp, op, dw, eflag, comp_0, comp_1, comp_2, fflag):    
        self.custom_energies = c_en
        self.icru_mat = ic_mat
        self.density_mat = den_mat
        self.mean_free_p = mfp
        self.op = op
        self.do_what = dw
        self.energy_flag = eflag
        self.frac_flag = fflag
        self.comp_0 = comp_0
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.weight_frac_list = []
        self.fetch_compound()
        self.calc_weight_fraction()
        
    def calc_weight_fraction(self):
        keys = self.dict_comp.keys()
        func = np.vectorize(lambda i : element(i).mass * self.dict_comp[i])
        keys = func(keys)
        self.weight_fraction = keys/np.sum(keys)
        if self.frac_flag:
            self.weight_fraction = np.asarray(self.weight_frac_list)

        if not self.frac_flag:
            values = self.dict_comp.values()
            self.number_fraction = values/np.sum(values)

    def comp_input(self):
        if self.frac_flag:
            za = self.comp_1
            weight_frac_list = self.comp_2.split()
            self.weight_frac_list[:] = [float(x) for x in weight_frac_list]
        else:
            za = self.comp_0
        return za
    
    def fetch_compound(self):
        """This function takes user input of the compound."""
        za = self.comp_input()
        compound_z_list = za.split()
        compound_z_list = [x.capitalize() for x in compound_z_list]
        za = ''.join(compound_z_list)
        sub = Substance.from_formula(za)
        self.sub_unicode_name = sub.unicode_name
        self.dict_comp = sub.composition
        self.ta_param = np.full((len(self.dict_comp),7,80), np.nan)
        self.formula_for_text = za

    def total_attenuation(self):
        """Calculates total attenuation coefficients from NIST database."""
        for i,z in enumerate(self.dict_comp.keys()):
            loc = 'NIST/'+'MDATX3n.' + str(z)
            m = element(z).mass
            data = np.loadtxt(loc, delimiter=', ').T
            data[1:6] * n * 10**(-24) / m
            data = np.append(data.T, [np.sum(data[1:6], axis = 0)], axis = 0)
            self.ta_param[i] = data
        return self.ta_param
    
    def d1(self):
        func = np.vectorize(lambda i : element(i).mass)
        keys = self.dict_comp.keys()
        denom = self.weight_fraction / func(keys) 
        denom1 = np.sum(denom) * n
        return denom1

    def myu(self):
        denom1 = self.d1()
        self.myu_comp = np.ndarray((7,80))
        params = self.total_attenuation()
        self.myu_comp = np.sum((params.T * self.weight_fraction).T, axis=0)
        self.myu_comp[0] = params[0][0]
        self.myu_comp = np.append(self.myu_comp, [self.myu_comp[-1]/denom1], axis = 0)
        self.myu_comp = np.append(self.myu_comp, [(self.myu_comp[-2]-self.myu_comp[1])/denom1], axis = 0)
        return self.myu_comp

    def zeff_by_log(self):
        photon_abs_element_list = loadtxt(
            "element_photo_abs", usecols=(0), unpack=True).reshape((99,80))
        
        if self.frac_flag:
            self.photon_comp = self.myu_comp[-2]
        else:        
            params = self.total_attenuation()
            func = np.vectorize(lambda i : element(i).mass)
            keys = self.dict_comp.keys()
            self.photon_comp = np.sum((params.T * self.number_fraction * func(keys)).T, axis=0)[-3] / n

        params = photon_abs_element_list.T    
        zno = np.arange(1,99)    
        z_comp = self.dict_comp.keys()
        zeff = np.full(80, np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        func = np.vectorize(lambda i : element(i).mass)
        Aavg = np.sum(self.dict_comp.values() * func(z_comp)) / np.sum(self.dict_comp.values())

        func = np.vectorize(lambda pa, ph : InterpolatedUnivariateSpline(zno, pa - ph).roots())
        func = np.vectorize(lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg)))
        zeff = func(params, self.photon_comp)

        return zeff

    def zeq_by_R(self, gp=False):
        R_element_list = loadtxt("element_R", usecols=(0), unpack=True).reshape((99,80))

        params = self.total_attenuation()
        R_comp_1 = np.sum((params.T * self.weight_fraction).T, axis=0)[2]
        R_comp_2 = np.sum((params.T * self.weight_fraction).T, axis=0)[-1]
        self.R_comp = R_comp_1/R_comp_2

        params = R_element_list.T    
        zno = np.arange(1,99)    
        z_comp = self.dict_comp.keys()
        self.zeq = np.full(80, np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        func = np.vectorize(lambda i : element(i).mass)
        Aavg = np.sum(self.dict_comp.values() * func(z_comp)) / np.sum(self.dict_comp.values())

        func = np.vectorize(lambda pa, ph : InterpolatedUnivariateSpline(zno, pa - ph).roots())
        func = np.vectorize(lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg)))
        self.zeq = func(params, self.R_comp)    
        #[np.nonzero(np.in1d(e1, e_range))]
        prob_flag = True
    
        if gp and prob_flag:
            b = get_gp('A',1)
            c = get_gp('A',2)
            a = get_gp('A',3)
            xk = get_gp('A',4)
            d = get_gp('A',5)

            b1 = get_gp('B',1)
            c1 = get_gp('B',2)
            a1 = get_gp('B',3)
            xk1 = get_gp('B',4)
            d1 = get_gp('B',5)
            
            B = {}
            BE = {}

            #mfp_q = prompt(user_mfp_list, style=style_1)['user_mfp_list']
            mfp_q = mean_free_path
            mfps = mfp_q.split()
            mfps[:] = [float(x) for x in mfps]
            for x in mfps:
                for i,e in enumerate(gp_energy_range):
                    k1 = (
                        (c[i] * (x**a[i])) + d[i] *( (np.tanh( (x/xk[i])-2 ) - np.tanh( -2 )) / ( 1 - np.tanh( -2 ) ))
                    )
                    if mfps.index(x) == 0:
                        B[e*1000000] = [round(b[i],3), round(c[i],3), round(a[i],3), round(xk[i],3), round(d[i],3) ]
                    if k1 == 1:
                        B[e*1000000].append(round(1 +(( b[i]-1 ) * x),3) )
                    else:
                        B[e*1000000].append(round(1 +( ( (b[i] - 1) * ( k1**x - 1) ) / ( k1 - 1) ),3) )


                for i,e in enumerate(gp_energy_range):
                    k1 = (
                        (c1[i] * (x**a1[i])) + d1[i] *( (np.tanh( (x/xk1[i])-2 ) - np.tanh( -2 )) / ( 1 - np.tanh( -2 ) ))
                    )
                    if mfps.index(x) == 0:
                        BE[e*1000000] = [ round(b1[i],3), round(c1[i],3), round(a1[i],3), round(xk1[i],3), round(d1[i],3) ]
                    if k1 == 1:
                        BE[e*1000000].append(round(1 +(( b1[i]-1 ) * x),3) )
                    else:
                        BE[e*1000000].append( round(1 + (( (b1[i] - 1) * ( k1**x - 1) ) / ( k1 - 1)),3) )
            return 'G-P fitting parameters and buildup factors - EABF, EBF'


    def get_gp(self, db, param):
        good_z = []
        all_y = np.full((22,25), np.nan)
        all_z = np.full(22, np.nan)
        
        for i,j in enumerate(range(4,83)):
            good_z.append(j)    
            zin = str(j)
            loc = 'ANSI_data/ANSI_'+db+'/DATAn_'+zin
            try:
                a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,param],dtype=float)
                y1 = array(a[param])
                x = array(a[0])
                all_y[i] = y1
                all_z[i] = j
            except FileNotFoundError:
                continue
            
        a = all_y.T
        b = self.zeq[np.nonzero(np.in1d(self.ta_param[0], x))]
        # new_good_z = np.append(b, all_z).sort()
        # a = list(zip(*all_y))
        # inter_y=[]
        # new_good_z = []
        # for e in gp_energy_range:
        #     b = list(all_z)
        #     b.append(zeq_R[e*1000000][0])
        #     b.sort()
        #     new_good_z.append(b)
        f = interp1d(all_z, a, kind = 'cubic')
        inter_y = f(b)
        return inter_y
        # for i in range(len(a)):
        #     f = interp1d(all_z, a[i], kind='cubic')
        #     inter_y = f(new_good_z[i])
        #     new_param.append(inter_y[new_good_z[i].index(zeq_R[gp_energy_range[i]*1000000][0])])
        # return new_param
        
def CreateFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

def CreateLog(dat1, loc = "InputLog.log"):
    with open(loc,"a") as text_file:
        text_file.write(json.dumps(dat1))
        text_file.write('\n')

CreateFolder('Excel_Sheets')


@eel.expose
def main1( comp_0a='H', do_what_now='Back', output='Do both', ff1=False,comp_1a='H', comp_2a='1', eflag=False,  mfp='1', density='1', \
    rel_mat='Air', custom_energies_list=''):
    comp = Compound(custom_energies_list, rel_mat, density, mfp, output, do_what_now, eflag, 
                    comp_0a, comp_1a, comp_2a, ff1)
    input_log = {}
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    input_log["Log time"] = dt_string
    input_log["Composition_0"] = comp_0
    input_log["Composition_1a"] = comp_1
    input_log["Composition_1b"] = comp_2
    input_log["Composition_1b"] = comp_2
    input_log["Parameter"] = do_what
    input_log["Output"] = op
    input_log["Energies"] = custom_energies
    input_log["KERMA Material"] = icru_mat
    input_log["Density"] = density_mat
    input_log["MFP"] = mean_free_path
    
    CreateLog(input_log)
    start_time = time.process_time()

    fetch_compound()
    myu()
    methods()
    CreateLog(f'Time elapsed: {time.process_time() - start_time}s')
    eel.excel_alert("Computation complete!")

eel.init('web')
eel.start('landing2.4.html',size=(1024, 550))