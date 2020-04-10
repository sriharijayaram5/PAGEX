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
        self.mean_free_path = mfp
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
        func = np.vectorize(lambda i : element(int(i)).mass * self.dict_comp[i])
        keys = func([*keys])
        self.weight_fraction = keys/np.sum(keys)
        if self.frac_flag:
            self.weight_fraction = np.asarray(self.weight_frac_list)

        if not self.frac_flag:
            values = np.asarray([*self.dict_comp.values()])
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
            data[1:] = data[1:] * n * 10**(-24) / m
            data = np.append(data, [np.sum(data[1:], axis = 0)], axis = 0)
            self.ta_param[i] = data
        return self.ta_param
    
    def d1(self):
        func = np.vectorize(lambda i : element(int(i)).mass)
        keys = self.dict_comp.keys()
        denom = self.weight_fraction / func([*keys]) 
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
        dest_filename = 'Save_File/Photon mass attenuation and interaction cross section parameters'
        if len(self.dict_comp)==1:
            header = ['Energy (MeV)', 'σ(coh)  (cm²/atom)','σ(incoh)  (cm²/atom)','σ(pe)  (cm²/atom)','σ(pair)  (cm²/atom)','σ(trip)  (cm²/atom)'
            ,'σ(wo/coh) (cm²/atom)','σ(w/coh) (cm²/atom)','μ/ρ(coh)  (cm²/g)','μ/ρ(incoh)  (cm²/g)','μ/ρ(pe)  (cm²/g)'
            ,'μ/ρ(pair)  (cm²/g)','μ/ρ(trip)  (cm²/g)','μ/ρ(wo/coh) (cm²/g)','μ/ρ(w/coh) (cm²/g)']
            func = np.vectorize(lambda i : element(int(i)).mass/n)
            keys = self.dict_comp.keys()
            f = func([*keys])
            params = [*(self.myu_comp[0:-3]*f), self.myu_comp[-1]*f, self.myu_comp[-3]*f, *self.myu_comp[1:-3], self.myu_comp[-1], self.myu_comp[-3]]
        else:
            header = ['Energy (MeV)', 'μ/ρ(coh)  (cm²/g)','μ/ρ(incoh)  (cm²/g)','μ/ρ(pe)  (cm²/g)','μ/ρ(pair)  (cm²/g)',
            'μ/ρ(trip)  (cm²/g)','μ/ρ(wo/coh) (cm²/g)','μ/ρ(w/coh) (cm²/g)']
            params = [*self.myu_comp[0:-3], self.myu_comp[-1], self.myu_comp[-3]]

        data = {'name' : dest_filename, 'header' : header,'params' : params}
        data['plot_params'] = [{'para_name' : '\dfrac{\mu}{\\rho}\ \ (cm^{2}/g)', 'value' : self.myu_comp[-3] }]
        return data

    def zeff_by_log(self):
        photon_abs_element_list = loadtxt(
            "element_photo_abs", usecols=(0), unpack=True).reshape((99,80))
        self.myu()
        func = np.vectorize(lambda i : element(int(i)).mass)
        keys = [*self.dict_comp.keys()]
        params = self.total_attenuation()
        self.photon_comp = np.sum((params.T * self.number_fraction * func(keys)).T, axis=0)[-1] / n
        if self.frac_flag:
            self.photon_comp = self.myu_comp[-2]
        params1 = photon_abs_element_list.T    
        zno = np.arange(1,100)    
        z_comp = [*self.dict_comp.keys()]
        zeff = np.full(80, np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(80):
            zeff[i] = func(params1[i], self.photon_comp[i])
        
        dest_filename = 'Save_File/Photon Zeff - Interpolation method'
        data = {'name' : dest_filename, 'header' : ['Energy (MeV)', 'Zeff', 'Neff (electrons/g)'], 'params' : [params[0][0], zeff, self.myu_comp[-3]/self.photon_comp*zeff]}
        data['plot_params'] = [{'para_name' : 'Z_{eff}', 'value' : zeff }]
        return data

    def zeq_by_R(self, gp=False):
        R_element_list = loadtxt("element_R", usecols=(0), unpack=True).reshape((99,80))

        params = self.total_attenuation()
        R_comp_1 = np.sum((params.T * self.weight_fraction).T, axis=0)[2]
        R_comp_2 = np.sum((params.T * self.weight_fraction).T, axis=0)[-1]
        self.R_comp = R_comp_1/R_comp_2

        params1 = R_element_list.T    
        zno = np.arange(1,100)    
        z_comp = [*self.dict_comp.keys()]
        self.zeq = np.full(80, np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(80):
            self.zeq[i] = func(params1[i], self.R_comp[i])    

        dest_filename = 'Save_File/Photon Zeq'
        data = {'name' : dest_filename, 'header' : ['Energy (MeV)', 'Zeq', 'R'], 'params' : [params[0][0], self.zeq, self.R_comp]}
        data['plot_params'] = [{'para_name' : 'Z_{eq}', 'value' : self.zeq }]
        if gp:
            b = self.get_gp('A',1)
            c = self.get_gp('A',2)
            a = self.get_gp('A',3)
            xk = self.get_gp('A',4)
            d = self.get_gp('A',5)

            b1 = self.get_gp('B',1)
            c1 = self.get_gp('B',2)
            a1 = self.get_gp('B',3)
            xk1 = self.get_gp('B',4)
            d1 = self.get_gp('B',5)
            
            self.p_B = [b, c, a, xk, d]
            self.p_BE = [b1, c1, a1, xk1, d1]
            
            mfp = [float(x) for x in self.mean_free_path.split()]
            self.B = np.full((len(mfp), len(b)), np.nan)
            self.BE = np.full((len(mfp), len(b)), np.nan)
            for l,mfps in enumerate(mfp):
                f = lambda x : (
                                    (c * (x**a)) + d *( (np.tanh( (x/xk)-2 ) - np.tanh( -2 )) / ( 1 - np.tanh( -2 ) ))
                                ) 
                k1 = np.asarray(f(mfps))
                f = lambda x : np.where(k1==1, 1 +(( b-1 ) * x), 1 +( ( (b - 1) * ( k1**x - 1) ) / ( k1 - 1) ))
                self.B[l] = f(mfps)
                f = lambda x : (
                                    (c1 * (x**a1)) + d1 *( (np.tanh( (x/xk1)-2 ) - np.tanh( -2 )) / ( 1 - np.tanh( -2 ) ))
                                ) 
                k1 = np.asarray(f(mfps))
                f = lambda x : np.where(k1==1, 1 +(( b1-1 ) * x), 1 +( ( (b1 - 1) * ( k1**x - 1) ) / ( k1 - 1) )) 
                self.BE[l] = f(mfps)
            
            dest_filename = 'Save_File/G-P fitting parameters and buildup factors'
            header = ['Energy (MeV)', 'b', 'c', 'a', 'Xk', 'd'] + [f'EABF {i}mfp' for i in mfp]
            header += ['b1', 'c1', 'a1', 'Xk1', 'd1'] + [f'EBF {i}mfp' for i in mfp]
            loc = 'ANSI_data/ANSI_'+'A'+'/DATAn_'+'4'
            a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0],dtype=float)
            energy = np.asarray(a[0])
            data = {'name' : dest_filename, 'header' : header, 'params' : [energy, *self.p_B, *self.B, *self.p_BE, *self.BE]}
            for i, m in enumerate(mfp):
                data['plot_params'] = [{'para_name' : f'EABF\ at\ {m}\ mfp', 'value' : self.B[i] },
                                        {'para_name' : f'EBF\ at\ {m}\ mfp', 'value' : self.BE[i] }]
        return data


    def get_gp(self, db, param):
        all_y = []
        all_z = []
        for j in range(4,83):
            zin = str(j)
            loc = 'ANSI_data/ANSI_'+db+'/DATAn_'+zin
            try:
                a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,param],dtype=float)
                y = a[param]
                x = a[0]
                all_y.append(y)
                all_z.append(j)
            except FileNotFoundError:
                continue
            
        a = np.asarray(all_y).T
        b = self.zeq[np.nonzero(np.in1d(self.ta_param[0][0], x*1000000))]
        inter_y = np.full_like(b, np.nan)
        for j,i in enumerate(a):    
            f = interp1d(all_z, i, kind = 'cubic')
            inter_y[j] = f(b[j])
        return inter_y
    
    def zeff_by_Ratio(self):
        '''Zeff by direct method'''
        params = self.total_attenuation()
        self.myu()
        func = np.vectorize(lambda i : element(int(i)).mass)
        keys = [*self.dict_comp.keys()]
        self.photon_comp = np.sum((params.T * self.number_fraction * func(keys)).T, axis=0)[-1] / n
        self.photon_e_comp = np.sum((params.T * self.number_fraction * func(keys) / keys).T, axis=0)[-1] / n
        self.zeff_ratio = self.photon_comp / self.photon_e_comp
    
        dest_filename = 'Save_File/Photon Zeff - Direct method'
        data = {'name' : dest_filename, 
                'header' : ['Energy (MeV)', 
                'σₐ Average Cross Section per Atom (cm²/atom)',
                'σₑ Average Cross Section per Electron (cm²/electron)','Zeff',
                'Neff (electrons/g)'], 'params' : [params[0][0],  self.photon_comp, self.photon_e_comp, self.zeff_ratio,self.myu_comp[-3]/self.photon_e_comp]}
        data['plot_params'] = [{'para_name' : 'Z_{eff}', 'value' : self.zeff_ratio }]
        return data

    def stopping_power_compound_post(self):      
        formula = []
        for i,e in enumerate(self.dict_comp.keys()):
            formula.extend(
                [element(e).symbol, str(round(self.weight_fraction[i],6)), '\n'])
        # formula=['H','0.11190674437968359','\n','O','0.8880932556203164']
        formula_req = ''
        for f in formula:
            formula_req = formula_req + f + ' '
        with requests.Session() as c:
            url = 'https://physics.nist.gov/cgi-bin/Star/estar-ut.pl'
            d = self.density_mat
            data1 = {'Name': 'Material',
                    'Density': d,
                    'Formulae': formula_req}
            a = True
            while a:
                try:
                    r = c.post(url, data=data1)
                    a = False
                except requests.exceptions.ConnectionError:
                    eel.excel_alert("Please check internet conncection. This parameter requires an active connection.")

            soup = BeautifulSoup(r.text, 'html.parser')
            results = soup.find('input').attrs['value']
        with requests.Session() as c:
            url = 'https://physics.nist.gov/cgi-bin/Star/e_table-ut.pl'
            pairnum = -1
            for f in formula:
                try:
                    float(f)
                    pairnum = pairnum + 1
                except ValueError:
                    pass
            lines = []
            for f in range(0, len(formula), 3):
                if formula[f] == '\n':
                    continue
                lines.append(formula[f] + ' ' + formula[f + 1])
            data1 = {'I': results,
                    'ShowDefault': 'on',
                    'Name': 'Material',
                    'Density': d,
                    'pairnum': pairnum}
            for l in lines:
                data1['line' + str(lines.index(l))] = l
            a = True
            while a:
                try:
                    r = c.post(url, data=data1)
                    a = False
                except requests.exceptions.ConnectionError:
                    eel.excel_alert("Please check internet conncection. This parameter requires an active connection.")
            soup = BeautifulSoup(r.text, 'html.parser')
            results = soup.find_all('pre')

            dat1 = results[0].contents[12:-1:2]
            dat1[:] = [x for x in dat1 if x != ' ']
            dat1[:] = [x for x in dat1 if x != '\x0c']
            dat1[:] = [x.strip() for x in dat1]
            a, b, = loadtxt(dat1, usecols=(0, 3), unpack=True)
            return a, b
            
    def zeff_electron_interaction(self):
        x, mass_stopping_power = self.stopping_power_compound_post()
        electron_int_cross = mass_stopping_power / self.d1()
        electron_msp = mass_stopping_power
        ele_electron_int_cross = np.full((98,81), np.nan)
        for i in range(1, 99):
            if i < 10:
                zin = '00' + str(i)
            else:
                zin = '0' + str(i)
            loc = 'EStar_data/DATA' + zin
            a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,3])
            y = a[3] / (n / element(i).mass)
            ele_electron_int_cross[i-1] = y

        params = ele_electron_int_cross.T    
        zno = np.arange(1,99)    
        z_comp = [*self.dict_comp.keys()]
        self.zeff_ele = np.full((len(x)), np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(len(x)):
            self.zeff_ele[i] = func(params[i], electron_int_cross[i]) 

        dest_filename = 'Save_File/Electron interaction parameters'
        data = {'name' : dest_filename, 'header' : 
                ['Energy (MeV)', 'S(E)/ρ (MeV cm²/g)', 
                 'Sc (MeV cm²/atom)', 'Zeff', 'Neff (electrons/g)'], 
                'params' : [x, electron_msp, electron_int_cross, self.zeff_ele, electron_msp / electron_int_cross * self.zeff_ele]}
        data['plot_params'] = [{'para_name' : 'S(E)/\\rho\ (MeV\ cm^{2}/g)', 'value' : electron_msp },
                                {'para_name' : 'Z_{eff}', 'value' : self.zeff_ele }]
        return data
        
    def zeff_proton_interaction(self):
        good_z = np.arange(1,93)
        all_y = []
        all_z = []
        for i in range(1, 93):
            if i < 10:
                zin = '00' + str(i)
            else:
                zin = '0' + str(i)
            loc = 'PStar_data/DATA' + zin
            try:
                a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,3])
                x = a[0]
                y = a[3]
            except OSError:
                continue
            all_y.append(y)
            all_z.append(i)

        a = np.asarray(all_y).T
        f = interp1d(all_z, a, kind = 'cubic')
        inter_y = f(good_z)
        all_new_y = inter_y.T

        msp = all_new_y[np.asarray([*self.dict_comp.keys()])-1]
        msp_comp = np.sum((msp.T * self.weight_fraction).T, axis=0)
        proton_int_cross = msp_comp / self.d1()
        
        func = np.vectorize(lambda i : element(int(i)).mass)
        b = n / func(np.arange(1,93))
        ele_proton_int_cross = (all_new_y.T / b).T

        params = ele_proton_int_cross.T    
        zno = np.arange(1,93)    
        z_comp = [*self.dict_comp.keys()]
        self.zeff_proton = np.full(len(x), np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(len(x)):
            self.zeff_proton[i] = func(params[i], proton_int_cross[i]) 

        dest_filename = 'Save_File/Proton interaction parameters'
        data = {'name' : dest_filename, 'header' : 
                ['Energy (MeV)', 'S(E)/ρ (MeV cm²/g)', 
                 'Sc (MeV cm²/atom)', 'Zeff', 'Neff (electrons/g)'], 
                'params' : [x, msp_comp, proton_int_cross, self.zeff_proton, msp_comp / proton_int_cross * self.zeff_proton]}
        data['plot_params'] = [{'para_name' : 'S(E)/\\rho\ (MeV\ cm^{2}/g)', 'value' : msp_comp },
                                {'para_name' : 'Z_{eff}', 'value' : self.zeff_proton }]
        return data
        
    def zeff_alpha_interaction(self):
        good_z = np.arange(1,93)
        all_y = []
        all_z = []
        for i in range(1, 93):
            if i < 10:
                zin = '00' + str(i)
            else:
                zin = '0' + str(i)
            loc = 'AStar_data/DATA' + zin
            try:
                a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,3])
                x = a[0]
                y = a[3]
            except OSError:
                continue
            all_y.append(y)
            all_z.append(i)

        a = np.asarray(all_y).T
        f = interp1d(all_z, a, kind = 'cubic')
        inter_y = f(good_z)
        all_new_y = inter_y.T

        msp = all_new_y[np.asarray([*self.dict_comp.keys()])-1]
        msp_comp = np.sum((msp.T * self.weight_fraction).T, axis=0)
        alpha_int_cross = msp_comp / self.d1()
        
        func = np.vectorize(lambda i : element(int(i)).mass)
        b = n / func(np.arange(1,93))
        ele_alpha_int_cross = (all_new_y.T / b).T

        params = ele_alpha_int_cross.T    
        zno = np.arange(1,93)    
        z_comp = [*self.dict_comp.keys()]
        self.zeff_alpha = np.full(len(x), np.nan)
        
        avg = np.sum(z_comp * self.weight_fraction)
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(len(x)):
            self.zeff_alpha[i] = func(params[i], alpha_int_cross[i]) 

        dest_filename = 'Save_File/Alpha particle interaction parameters'
        data = {'name' : dest_filename, 'header' : 
                ['Energy (MeV)', 'S(E)/ρ (MeV cm²/g)', 
                 'Sc (MeV cm²/atom)', 'Zeff', 'Neff (electrons/g)'], 
                'params' : [x, msp_comp, alpha_int_cross, self.zeff_alpha, msp_comp / alpha_int_cross * self.zeff_alpha]}
        data['plot_params'] = [{'para_name' : 'S(E)/\\rho\ (MeV\ cm^{2}/g)', 'value' : msp_comp },
                                {'para_name' : 'Z_{eff}', 'value' : self.zeff_alpha }]
        return data

    def kerma_1(self, relative_to_choice='AIR', kerma = False):
        ele_x = []
        for i in range(1, 93):
            z = element(i).atomic_number
            zin = ''
            if z < 10:
                zin = '0' + str(z)
            else:
                zin = str(z)
            loc = 'XRay_data1/DATAn' + zin
            a = read_csv(loc,delim_whitespace=True,header=None,usecols=[0,2])
            x = a[0]
            y = a[2]
            ele_x.append(list(y))
        ele_x = np.asarray(ele_x)

        mec_comp = ele_x[np.asarray([*self.dict_comp.keys()])-1]
        self.mec_comp1 = np.sum((mec_comp.T * self.weight_fraction).T, axis=0)
        sigmaa = self.mec_comp1 / self.d1()

        loc = 'XRay_Comp1/DATAn_' + relative_to_choice
        self.mec_rel = loadtxt(loc, usecols=(2), unpack=True)
        
        self.kerma = self.mec_comp1 / self.mec_rel

        ele_x_int_cross = ele_x / self.d1()
        params = ele_x_int_cross.T    
        zno = np.arange(1,93)    
        z_comp = [*self.dict_comp.keys()]
        self.zeff_x = np.full(len(x), np.nan)

        avg = np.sum(z_comp * self.weight_fraction)
        func = lambda pa, ph : min(InterpolatedUnivariateSpline(zno, pa - ph).roots(), key = lambda x : abs(x - avg))
        for i in range(len(x)):
            self.zeff_x[i] = func(params[i], sigmaa[i]) 
    
        if kerma:
            dest_filename = 'Save_File/Relative KERMA'
            data = {'name' : dest_filename, 'header' : 
                    ['Energy (MeV)', f'{relative_to_choice} KERMA'], 
                    'params' : [x, self.kerma]}
            data['plot_params'] = [{'para_name' : f'{relative_to_choice}\ KERMA', 'value' : self.kerma }]
        else:    
            dest_filename = 'Save_File/Photon mass-energy absorption coefficients'
            data = {'name' : dest_filename, 'header' : 
                    ['Energy (MeV)', 'MEC μₑₙ/ρ (cm²/g)', 
                    'Z PEAeff', 'N PEAeff (electrons/g)'], 
                    'params' : [x, self.mec_comp1, self.zeff_x, self.mec_comp1/sigmaa*self.zeff_x]}
            data['plot_params'] = [{'para_name' : '\mu_{en}/\\rho\ (cm^{2}/g)', 'value' : self.mec_comp1 }]
        return data
     
    def write_to_csv(self, data):
        fname = data['name'] + f'-{self.formula_for_text}.csv'
        X = np.asarray(data['params'])    
        if X[0][0] > 1:
                X[0] = X[0]/10e6      
        np.savetxt(fname, X.T, header = ', '.join(data['header']), delimiter = ', ', comments = '#', fmt = '%s', encoding = "utf-8")
    
    def plot_parameter(self, data):
        x = data['params'][0]
        plot_params = data['plot_params']
        for para in plot_params:    
            plt.ylabel('$%s$'%para['para_name'], fontname='Calibri')
            plt.xlabel('$E\ (MeV)$', fontname='Calibri')
            plt.ticklabel_format(
                axis='both', style='sci', scilimits=(
                    0,0), useMathText=True)
            plt.tick_params(
                axis='both',direction='in',which='both',top=True,right=True
            )
            if x[0] > 1:
                x = x/10e6
            name = Substance.from_formula(self.formula_for_text).unicode_name
            if self.frac_flag:
                    sub = name +'('
                    for j in self.weight_fraction:
                        sub = sub + str(round(j,2)) + ','
                    sub = sub.strip(',')
                    name = sub + ')'
            if para['para_name'] in ['Z_{eff}','Z_{eq}','Relative\KERMA','Z_{PEAeff}']:
                plt.semilogx(x, para['value'], '-x', markersize=5, label=name)
            else:
                plt.loglog(x, para['value'], '-x', markersize=5, label=name)
            plt.legend(loc='upper right')
            plt.show()
            plt.close()

    def interpolate_e(self, energy, parameter, new_energy, num=0):
        y=[]
        
        if isinstance(parameter[energy[0]],float): 
            for j in energy:
                y.append(parameter[j])
        else:
            for j in energy: 
                y.append(parameter[j][num])
        f = interp1d(energy,y,kind='cubic')
        return(f(new_energy))

    def interpolate_e_linear(self, energy, parameter, new_energy, num=0):
        y=[]
        
        if isinstance(parameter[energy[0]],float): 
            for j in energy:
                y.append(parameter[j])
        else:
            for j in energy: 
                y.append(parameter[j][num])
        f = interp1d(energy,y,kind='slinear')
        return(f(new_energy))

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

CreateFolder('Save_File')


@eel.expose
def main1( comp_0a=None, do_what_now=None, output=None, ff1=False, comp_1a=None, comp_2a=None, eflag=False,  mfp=None, density=None, \
    rel_mat=None, custom_energies_list=None):
    comp = Compound(custom_energies_list, rel_mat, density, mfp, output, do_what_now, eflag, 
                    comp_0a, comp_1a, comp_2a, ff1)
    input_log = {}
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    input_log["Log time"] = dt_string
    if ff1:
        input_log["Elements"] = comp_1a
        input_log["Weight Fraction"] = comp_2a
    else:
        input_log["Compound"] = comp_0a
    input_log["Parameter"] = do_what_now
    if not eflag:
        input_log["Energies"] = custom_energies_list
    if rel_mat is not None:
        input_log["KERMA Material"] = rel_mat
    if density is not None:
        input_log["Density"] = density
    if mfp is not None:
        input_log["MFP"] = mfp
    param = do_what_now
    params = [
        'Partial/total interaction cross sections and mass attenuation coefficients',
        'Photon energy absorption coefficients (cm²/g), Z PEAeff, N PEAeff (electrons/g)',
        'Relative KERMA',
        'Equivalent atomic number - Zeq',
        'G-P fitting parameters and buildup factors - EABF, EBF',
        'Direct method',
        'Interpolation method',
        'Proton interaction',
        'Electron interaction',
        'Alpha particle interaction'
    ]
    output_choices = [

        'Write parameter to excel sheet only',
        'Plot Energy vs. Parameter only',
        'Do both'
    ]
    CreateLog(input_log)
    start_time = time.process_time()
    if param == params[0]:
        data = comp.myu()
    elif param == params[1]:
        data = comp.kerma_1()
    elif param == params[2]:
        data = comp.kerma_1(relative_to_choice =  rel_mat, kerma = True)
    elif param == params[3]:
        data = comp.zeq_by_R()
    elif param == params[4]:
        data = comp.zeq_by_R(gp = True)
    elif param == params[5]:
        data = comp.zeff_by_Ratio()
    elif param == params[6]:
        data = comp.zeff_by_log()
    elif param == params[7]:
        data = comp.zeff_proton_interaction()
    elif param == params[8]:
        data = comp.zeff_electron_interaction()
    elif param == params[9]:
        data = comp.zeff_alpha_interaction()
    CreateLog(f'Time elapsed: {time.process_time() - start_time}s')
    if output == output_choices[0]:
        comp.write_to_csv(data)
    elif output == output_choices[1]:
        comp.plot_parameter(data)
    else:
        comp.write_to_csv(data)
        comp.plot_parameter(data)
    del(comp)
    # eel.excel_alert("Computation complete!")

eel.init('web')
eel.start('landing2.4.html',size=(1024, 550))