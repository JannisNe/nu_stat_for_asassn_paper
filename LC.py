import re # use regular patters
import sys # system commands
import string as string # string functions4
import math
import numpy as np # numerical tools
from scipy import *
import os, io
import pandas as pd
from matplotlib.pyplot import figure, show, rcParams
import matplotlib.pyplot as plt

# plt.rcParams['legend.numpoints']=1
# plt.rcParams['xtick.major.size'] = 11
# plt.rcParams['xtick.minor.size'] = 5
# plt.rcParams['ytick.major.size'] = 11
# plt.rcParams['ytick.minor.size'] = 5
# plt.rcParams['text.latex.preamble'] = [r'\boldmath']
# params = {'text.usetex': True, 'mathtext.fontset': 'stixsans'}
# rcParams.update(params)

#event
names=np.array(['AT2020rng','AT2020idu','AT2019fxr','AT2019dsg','SN2019aah','AT2018dyd','SN2018coq','AT2017hzv'])
nu_names = np.array(["IC210608A", "IC200530A", "IC200410A", "IC191001A", "IC191119A", "IC191001A", "IC190619A", "IC180908A"])
types = np.array(['Unknown', 'Dwarf Nova', 'Unknown', 'TDE', 'SN II', 'Unknown', 'SN II', 'Unknown'])
ra=np.array([333.1687755,253.335125,242.066735231,314.2623926,225.377973463,314.0556855,343.776167,144.801206241])
dec=np.array([+18.0510409,+27.435711,+6.82921933663,14.2044063,+4.87828907113,+12.8536367,+11.257611,-1.3396342480])
JDasassn=np.array([2459089.9,2458962.9,2458634.9,2458618.9,2458519.6,2458316.8,2458286.1,2458070.8])
JDicecube=np.array([2459373.7,2458999.8,2458950.5,2458758.3,2458806.5,2458758.3,2458654.1,2458370.3])

candidates_lc_dir = '../figures/candidates_lightcurves'

corrupted = ['coadd_bq287843', 'coadd_bq286657', 'coadd_bE242641']

def get_lc(name):
    
    m = names == name
    jdass = JDasassn[m]
    jdice = JDicecube[m]
    nu_name = nu_names[m]
    itype = types[m]
    data_dict = dict()
    
    for band in ['g', 'V']:
        with open(os.path.join('data/LC', f'{name}_{band}.dat'), 'r') as f:
            buf = io.StringIO(f.read().split('### ')[-1].replace('flux (mJy)', 'flux(mJy)'))
        try:
            data = pd.read_csv(buf, delim_whitespace=True, comment='#')
            data['mag'] = np.array([float(m.split('>')[-1]) for m in data.mag])
            corr_m = np.array([i in corrupted for i in data.IMAGE])
    #         m = data['flux(mJy)'] / data['flux_err'] >= sigma
    #         data = data[m]
            data_dict[band] = data[~corr_m]
            if np.any(corr_m):
                data_dict[f"corrupted_{band}"] = data[corr_m]
            
        except pd.errors.EmptyDataError:
            print(f"no data for {name} for {band}")
        
    return jdass, jdice, data_dict, nu_name[0], itype[0]
    
#     data=open('data/LC/%s_g.dat'%name,'r')
#     lines=data.readlines()
#     JD_g=[]
#     mag_g=[]
#     emag_g=[]
#     for j in range(size(lines)):
#         if (lines[j].split()[0]!='###'):

#             JD_g.append(float(lines[j].split()[0])-jdass)
#             if lines[j].split()[7][0]=='>':
#                 mag_g.append(float(lines[j].split()[7][1:]))
#             else:
#                 mag_g.append(float(lines[j].split()[7]))
#             emag_g.append(float(lines[j].split()[8]))


#     data=open('data/LC/%s_V.dat'%name,'r')
#     lines=data.readlines()
#     JD_V=[]
#     mag_V=[]
#     emag_V=[]
#     for j in range(size(lines)):
#         if (lines[j].split()[0]!='###'):

#             JD_V.append(float(lines[j].split()[0])-jdass)
#             if lines[j].split()[7][0]=='>':
#                 mag_V.append(float(lines[j].split()[7][1:]))
#             else:
#                 mag_V.append(float(lines[j].split()[7]))
#             emag_V.append(float(lines[j].split()[8]))
            
#     return jdass, jdice, JD_g, mag_g, emag_g, JD_V, mag_V, emag_V, nu_name[0], itype[0]


if __name__ == '__main__':
    
    # not tested!

    for name in names:

        jdass, jdice, JD_g, mag_g, emag_g, JD_V, mag_V, emag_V = get_lc(name)

        fig, ax1 = plt.subplots(figsize=(8,8), facecolor='w', edgecolor='k')  # ax1 is apparent magnitude
        plt.gca().invert_yaxis()
        plt.minorticks_on()

        #Photo
        #non detection
        ax1.plot(np.array(JD_V)[np.array(emag_V)==99.990],np.array(mag_V)[np.array(emag_V)==99.990],'^g', ls='none',alpha=0.5)

        ax1.errorbar(np.array(JD_V)[np.array(emag_V)!=99.990],np.array(mag_V)[np.array(emag_V)!=99.990],yerr=np.array(emag_V)[np.array(emag_V)!=99.990],color='g',marker='s', ls='none',alpha=0.8,label=r'\textbf{V}')

        ax1.plot(np.array(JD_g)[np.array(emag_g)==99.990],np.array(mag_g)[np.array(emag_g)==99.990],'^b', ls='none',alpha=0.5)
        ax1.errorbar(np.array(JD_g)[np.array(emag_g)!=99.990],np.array(mag_g)[np.array(emag_g)!=99.990],yerr=np.array(emag_g)[np.array(emag_g)!=99.990],color='b',marker='s', ls='none',alpha=0.8,label=r'\textbf{g}')

        ax1.axvline(x=jdice-jdass,color='r',label=r'\textbf{IceCube alert}')

        ax1.tick_params(axis='both', which='major', labelsize=25)
        ax1.set_xlabel(r'\textbf{JD - %s [days]}'%jdass,fontsize=25,fontweight='bold')
        ax1.set_ylabel(r'\textbf{Apparent magnitude [mag]}',fontsize=25,fontweight='bold')
        #ax1.set_ylim([21,17])

        ax1.legend(loc=0,title='',markerscale=1.5,ncol=1,prop={'size':18})
        ax1.minorticks_on()
        plt.title(r'\textbf{%s}'%(name), fontsize=20)
        plt.savefig('Figures_LC/%s.png'%name,bbox_inches='tight')
        plt.savefig('Figures_LC/%s.pdf'%name,dpi=600,bbox_inches='tight')
        plt.show()
        plt.close("all")
