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
    
    def myu(self):
        func = np.vectorize(lambda i : element(i).mass)
        keys = self.dict_comp.keys()
        denom = self.weight_fraction / func(keys) 
        denom1 = np.sum(denom) * n
        self.myu_comp = np.ndarray((7,80))
        params = self.total_attenuation()
        self.myu_comp = np.sum((params.T * self.weight_fraction).T, axis=0)
        self.myu_comp[0] = params[0]
        self.myu_comp = np.append(self.myu_comp, [self.myu_comp[-1]/denom1], axis = 0)
        self.myu_comp = np.append(self.myu_comp, [(self.myu_comp[-2]-self.myu_comp[1])/denom1], axis = 0)
        return self.myu_comp

        
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