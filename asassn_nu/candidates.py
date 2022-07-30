import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerTuple
from matplotlib.container import Container
import pandas as pd
import numpy as np
import os, requests, tarfile, zipfile, io
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u, io
# from style import output_folder, big_fontsize, base_width, base_height, dpi
import seaborn as sns
import json
from astropy.time import Time

from asassn_nu.info import output_folder_figures
from style import base_width, base_height, dpi, bandmark, bandcols, lc_uplim_kw, lc_errbar_kw, big_fontsize


candidates_lc_dir = os.path.join(output_folder_figures, '/candidates_lightcurves')

#########################################################
#           Candidate Counterparts
#########################################################

names = np.array(
    ['AT2020rng', 'AT2020idu', 'AT2019fxr', 'AT2019dsg', 'SN2019aah', 'AT2018dyd', 'SN2018coq', 'AT2017hzv'])
nu_names = np.array(
    ["IC210608A", "IC200530A", "IC200410A", "IC191001A", "IC191119A", "IC191001A", "IC190619A", "IC180908A"])
types = np.array(['Unknown', 'Dwarf Nova', 'Unknown', 'TDE', 'SN II', 'Unknown', 'SN II', 'Unknown'])
ra = np.array(
    [333.1687755, 253.335125, 242.066735231, 314.2623926, 225.377973463, 314.0556855, 343.776167, 144.801206241])
dec = np.array(
    [+18.0510409, +27.435711, +6.82921933663, 14.2044063, +4.87828907113, +12.8536367, +11.257611, -1.3396342480])
JDasassn = np.array([2459089.9, 2458962.9, 2458634.9, 2458618.9, 2458519.6, 2458316.8, 2458286.1, 2458070.8])
JDicecube = np.array([2459373.7, 2458999.8, 2458950.5, 2458758.3, 2458806.5, 2458758.3, 2458654.1, 2458370.3])

#########################################################
#           Bad images
#########################################################

corrupted = ['coadd_bq287843', 'coadd_bq286657', 'coadd_bE242641']

#########################################################
#           Get lightcurve data
#########################################################


def get_lc(name):
    m = names == name
    jdass = JDasassn[m]
    jdice = JDicecube[m]
    nu_name = nu_names[m]
    itype = types[m]
    data_dict = dict()

    for band in ['g', 'V']:
        with open(os.path.join('../data/LC', f'{name}_{band}.dat'), 'r') as f:
            buf = io.StringIO(f.read().split('### ')[-1].replace('flux (mJy)', 'flux(mJy)'))
        try:
            data = pd.read_csv(buf, delim_whitespace=True, comment='#')
            data['mag'] = np.array([float(m.split('>')[-1]) for m in data.mag])
            corr_m = np.array([i in corrupted for i in data.IMAGE])
            data_dict[band] = data[~corr_m]
            if np.any(corr_m):
                data_dict[f"corrupted_{band}"] = data[corr_m]

        except pd.errors.EmptyDataError:
            print(f"no data for {name} for {band}")

    return jdass, jdice, data_dict, nu_name[0], itype[0]


#########################################################
#           Plot lightcurve
#########################################################


def plot_lc(name, sigma=5):
    jdass, jdice, data, nu_name, ttype = get_lc(name)
    fig, ax = plt.subplots(figsize=(base_width, base_height), dpi=dpi)
    conts = list()
    labels = list()
    for band, d in data.items():

        if 'corrupted' in band:
            continue

        jd = d['JD'] - jdass
        mag = d['mag']
        err = d['mag_err']

        if len(jd):
            errm = (d['flux(mJy)'] < d.flux_err * sigma) | (err == 99.99)
            cont1 = ax.errorbar(jd[~errm], mag[~errm], yerr=err[~errm], color=bandcols[band], marker=bandmark[band],
                                **lc_errbar_kw)
            cont2, = ax.plot(jd[errm], mag[errm], color=bandcols[band], markeredgecolor=bandcols[band], **lc_uplim_kw)
            conts.append((cont1, cont2))
            labels.append(f"{band}, upper limits")

    ylim = ax.get_ylim()
    ax.set_ylim([ylim[1], ylim[0]])
    ax.set_ylabel('magnitude')
    mjd_ass = Time(jdass[0], format='jd').mjd
    mjd_ice = Time(jdice[0], format='jd').mjd
    ax.set_xlabel(f'MJD - {mjd_ass:.0f}')
    line = ax.axvline(mjd_ice - mjd_ass)
    ax.legend([line] + conts, [nu_name] + labels, handler_map={tuple: HandlerTuple(ndivide=None)})
    ax.set_title(name)
    fig.tight_layout()
    return fig, ax


def make_candidates_plot():
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(base_width * 2, base_height * 3), dpi=dpi, sharey='all',
                            gridspec_kw={'hspace': 0.35})
    sigma = 5
    ylims = list()
    inds = [0, 2, 3, 4, 6, 7]
    for i, j in enumerate(inds):
        name = names[j]
        ax = axs.flatten()[i]
        if name == 'AT2020idu':
            continue

        print(name)
        jdass, jdice, data, nu_name, ttype = get_lc(name)
        conts = list()
        labels = list()
        for band, d in data.items():
            if 'corrupted' in band:
                continue

            jd = d['JD'] - jdass
            mag = d['mag']
            err = d['mag_err']

            if len(jd):
                errm = (d['flux(mJy)'] < d.flux_err * sigma) | (err == 99.99)
                cont1 = ax.errorbar(jd[~errm], mag[~errm], yerr=err[~errm], color=bandcols[band], marker=bandmark[band],
                                    **lc_errbar_kw)
                cont2, = ax.plot(jd[errm], mag[errm], color=bandcols[band], **lc_uplim_kw)
                conts.append((cont1, cont2))
                labels.append(f"{band}, upper limits")

        if name == '_AT2020rng':
            t, mag, mage, fltr = get_ztffps_at2020rng()
            ulm = mage == 99.99
            for f in ['r', 'g']:
                m = fltr == f'ZTF_{f}'
                m = m & (t > -15) & (t < -11)
                kw = dict(lc_errbar_kw)
                kw['ecolor'] = None
                kw['markeredgecolor'] = None
                if np.any(m & ~ ulm):
                    ax.errorbar(t[m & ~ ulm], mag[m & ~ ulm], yerr=mage[m & ~ ulm], color=bandcols[f],
                                marker=bandmark[f], label=f"ZTF {f}", **kw)
                ax.scatter(t[m & ulm], mag[m & ulm], color=bandcols[f], marker='d', s=2, alpha=1,
                           label=f"ZTF {f} upper limits")

            disc_time = Time('2020-08-15 08:06:00.000').jd - jdass
            ax.scatter(disc_time, 16.6717, marker='o', color=bandcols['r'], facecolors='none', s=70,
                       label='ZTF discovery')
            ax.legend(ncol=2, loc='lower center')

        ylims.append(ax.get_ylim())
        mjd_ass = Time(jdass[0], format='jd').mjd
        mjd_ice = Time(jdice[0], format='jd').mjd
        line = ax.axvline(mjd_ice - mjd_ass)
        ax.text(.5, 1.03, f'{name} ({ttype})',
                horizontalalignment='center',
                transform=ax.transAxes, fontsize=big_fontsize - 2)
    ylims = np.array(ylims)
    axs[0][0].set_ylim(max(ylims[:, 1]), min(ylims[:, 0]))

    fig.supxlabel('Days after discovery')
    fig.supylabel('Magnitude')

    fig.legend([line] + conts, ['neutrino'] + labels, handler_map={tuple: HandlerTuple(ndivide=None)},
               bbox_to_anchor=(0.5, 0.95),
               bbox_transform=fig.transFigure, loc='upper center', ncol=3)

    fn = f'{candidates_lc_dir}/combined.pdf'
    print("saving under", fn)
    fig.savefig(fn)

    plt.show()
    plt.close()


#########################################################
#           AT2020rng
#########################################################

def get_ztffps_at2020rng():
    tab = io.ascii.read('data/at2020rng_ztffps.txt')
    names = list(tab[0])
    names[6] = 'info'
    names = [i.strip(',') for i in names]
    tab = io.ascii.read('data/at2020rng_ztffps.txt', names=names)[1:]
    tab = tab[tab['programid'] == '1']
    tab = tab[~(tab['forcediffimflux'] == 'null')]

    jdass, jdice, data, nu_name, ttype = get_lc('AT2020rng')
    f = np.array(tab['forcediffimflux']).astype(float)
    fu = np.array(tab['forcediffimfluxunc']).astype(float)
    zp = np.array(tab['zpdiff']).astype(float)
    t = np.array(tab['jd']).astype(float) - jdass
    fltr = tab['filter']

    m = f / fu > 5
    print(f"{sum(m)} detections")
    tab = tab

    mag = np.zeros(len(tab))
    mage = np.zeros(len(tab))

    mag[m] = zp[m] - 2.5 * np.log10(f[m])
    mage[m] = 1.0857 * fu[m] / f[m]

    mag[~m] = zp[~m] - 2.5 * np.log10(5 * fu[~m])
    mage[~m] = 99.99
    return t, mag, mage, fltr


def make_at2020rng_plot():

    jdass, jdice, data, nu_name, ttype = get_lc('AT2020rng')
    fig, ax = plt.subplots(figsize=(base_width * 2, base_height), dpi=dpi)
    conts = list()
    labels = list()
    sigma = 5
    for band, d in data.items():

        if 'corrupted' in band:
            continue

        jd = d['JD'] - jdass
        mag = d['mag']
        err = d['mag_err']

        if len(jd):
            errm = (d['flux(mJy)'] < d.flux_err * sigma) | (err == 99.99)
            cont1 = ax.errorbar(jd[~errm], mag[~errm], yerr=err[~errm], color=bandcols[band], marker=bandmark[band],
                                **lc_errbar_kw)
            cont2, = ax.plot(jd[errm], mag[errm], color=bandcols[band], markeredgecolor=bandcols[band], **lc_uplim_kw)
            conts.append((cont1, cont2))
            labels.append(f"ASAS-SN {band}, upper limits")

    t, mag, mage, fltr = get_ztffps_at2020rng()
    ulm = mage == 99.99

    for f in ['r', 'g']:
        m = fltr == f'ZTF_{f}'
        m = m & (t > -15) & (t < 2000)
        kw = dict(lc_errbar_kw)
        kw['ecolor'] = 'k'
        kw['markeredgecolor'] = 'k'
        if np.any(m & ~ ulm):
            marker = {'g': 'o', 'r': 'D'}[f]
            cont1 = ax.errorbar(t[m & ~ ulm], mag[m & ~ ulm], yerr=mage[m & ~ ulm], color=bandcols[f], marker=marker,
                                label=f"ZTF {f}", **kw)
        else:
            cont1 = None
        ztf_ulim_kw = dict(lc_uplim_kw)
        ztf_ulim_kw['marker'] = 'd'
        cont2, = ax.plot(t[m & ulm], mag[m & ulm], color=bandcols[f], label=f"ZTF {f} upper limits",
                         **ztf_ulim_kw)
        if cont1:
            conts.append((cont1, cont2))
            labels.append(f"ZTF {f}, upper limits")
        else:
            conts.append((cont2))
            labels.append(f"ZTF {f} upper limits")

    disc_time = Time('2020-08-15 08:06:00.000').jd - jdass
    cont4, = ax.plot(disc_time, 16.77, marker='o', color=bandcols['r'], mfc='none', ms=8, label='ZTF discovery', ls='')
    conts.append(cont4)
    labels.append(f"ZTF discovery epoch")

    ylim = ax.get_ylim()
    ax.set_ylim([ylim[1], ylim[0]])
    ax.set_ylabel('magnitude')
    mjd_ass = Time(jdass[0], format='jd').mjd
    mjd_ice = Time(jdice[0], format='jd').mjd
    ax.set_xlabel(f'MJD - {mjd_ass:.0f}')
    line = ax.axvline(mjd_ice - mjd_ass)
    ax.legend([line] + conts, [nu_name] + labels, handler_map={tuple: HandlerTuple(ndivide=None)},
              bbox_to_anchor=(0.5, 0.99), bbox_transform=fig.transFigure, loc='upper center', ncol=5)
    ax.set_title('')
    fig.tight_layout()

    ax.set_ylim([22, 16])
    plt.show()

    fn = f"{candidates_lc_dir}/at2020rng.pdf"
    print("saving under", fn)
    fig.savefig(fn, bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    make_candidates_plot()
    make_at2020rng_plot()
