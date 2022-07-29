import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerTuple
import pandas as pd
import numpy as np
import os
import io
from astropy.time import Time

from asassn_nu.info import data_dir, output_folder_figures
from asassn_nu.style import (
    base_width,
    base_height,
    dpi,
    lc_uplim_kw,
    bandmark,
    lc_errbar_kw,
    bandcols
)
from asassn_nu.fermi_coincidences import get_fermi_matches


fermi_lightcurves_dir = os.path.join(data_dir, 'fermi_matches_lightcurves')
fermi_lightcurves_plot_dir = os.path.join(output_folder_figures, 'fermi_lightcurves')


def plot_lc(r, unit='flux', sigma=5, timerange=None):
    if isinstance(timerange, type(None)):
        timerange = [-np.inf, np.inf]

    fgl_name = r['4FGL Name']
    counterpart_name = r['Counterpart']
    neutrino_name = r.name[0]
    neutrino_time_pre = r.name[0].strip('A').strip('B').split('IC')[-1]
    neutrino_mjd = Time(f"20{neutrino_time_pre[:2]}-{neutrino_time_pre[2:4]}-{neutrino_time_pre[-2:]}").mjd
    timerange = np.array(timerange) + neutrino_mjd

    fig, ax = plt.subplots(figsize=(base_width, base_height), dpi=dpi)
    if unit == 'flux':
        ax.axhline(0, ls='-', color='grey', alpha=1, lw=1)
    label_tot = list()
    conts_tot = list()

    for band in ['g', 'V']:

        with open(os.path.join(fermi_lightcurves_dir, f'{fgl_name.strip(" ").replace(" ", "_")}_{band}.dat'), 'r') as f:
            buf = io.StringIO(f.read().split('### ')[-1].replace('flux (mJy)', 'flux(mJy)'))
        data = pd.read_csv(buf, delim_whitespace=True, comment='#')
        mjdpre = np.array([Time(jd, format='jd').mjd for jd in data.JD])
        time_mask = (mjdpre > min(timerange)) & (mjdpre < max(timerange))
        data = data[time_mask]
        mjd = np.array([Time(jd, format='jd').mjd for jd in data.JD])

        alpha = np.zeros(len(data))

        if unit == 'flux':
            masks = [data['flux(mJy)'] < 5 * data.flux_err,
                     data['flux(mJy)'] >= 5 * data.flux_err]
            markers = [lc_uplim_kw['marker'], bandmark[band]]
            conts = list()
            lws = [dict(lc_uplim_kw), dict(lc_errbar_kw)]

            for mask, marker, kw in zip(masks, markers, lws):
                if not np.any(mask):
                    continue

                idata = data[mask]
                imjd = mjd[mask]
                kw['marker'] = marker
                cont = ax.errorbar(imjd, idata['flux(mJy)'], yerr=idata.flux_err, c=bandcols[band], **kw)
                conts.append(cont)

            label_tot.append(f"{band}, subthreshold")
            conts_tot.append((conts[-1], conts[0]))

        if unit == 'mag':

            uplim_mask = (data['flux(mJy)'] < data.flux_err * sigma) | (data['mag_err'] == 99.99)
            masks = [uplim_mask, ~uplim_mask]
            conts = list()

            mags = np.array([float(m.split('>')[-1]) for m in data.mag])

            l = ""
            if np.any(masks[1]):
                conts.append(ax.errorbar(mjd[masks[1]], mags[masks[1]], yerr=data.mag_err[masks[1]],
                                         marker=bandmark[band], color=bandcols[band], **lc_errbar_kw))
                l = band
            if np.any(masks[0]):
                conts.append(ax.plot(mjd[masks[0]], mags[masks[0]], color=bandcols[band], **lc_uplim_kw)[0])
                if l:
                    l = f"{band}, upper limits"
                else:
                    l = f'{band} upper limits'

            if np.any(masks[1]) or np.any(masks[0]):
                conts_tot.append(tuple(conts))
                label_tot.append(l)

    title = fgl_name if not isinstance(counterpart_name, str) else counterpart_name

    if unit == 'mag':
        ylim = ax.get_ylim()
        ax.set_ylim((ylim[-1], ylim[0]))

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_title(title.strip(' '), loc='center')
    line = ax.axvline(neutrino_mjd)
    ax.legend([line] + conts_tot, [neutrino_name] + label_tot, handler_map={tuple: HandlerTuple(ndivide=None)})
    ax.set_xlabel('MJD')
    xlabel = 'Flux [mJy]' if unit == 'Flux' else 'Magnitude'
    ax.set_ylabel(xlabel)

    return fig, ax, neutrino_mjd


def plot_fermi_matches_lightcurves():

    #########################################################
    #           Loading Fermi Matches
    #########################################################

    fermi_matches = get_fermi_matches()

    for i, r in fermi_matches.iterrows():
        fgl_name = r['4FGL Name']
        fn = [ifn for ifn in os.listdir(fermi_lightcurves_dir) if
              ifn.startswith(fgl_name.replace(' ', '_')) and ifn.endswith('.dat')]
        print(fgl_name, len(fn), fn)

        #########################################################
        #           Plot fluxes
        #########################################################

        fig, ax, nmjd = plot_lc(r)
        fn = os.path.join(fermi_lightcurves_plot_dir, f"{r.name[0]}_{r['4FGL Name'].replace(' ', '_')}.pdf")
        print(f"saving under {fn}")
        fig.subplots_adjust(bottom=0.2)
        fig.subplots_adjust(left=0.15)
        fig.savefig(fn)

        ax.set_xlim([nmjd - 100, nmjd + 100])
        fn = os.path.join(fermi_lightcurves_plot_dir, f"{r.name[0]}_{r['4FGL Name'].replace(' ', '_')}_closeup.pdf")
        print(f"saving under {fn}")
        fig.savefig(fn)
        plt.close()

        #########################################################
        #           Plot magnitudes
        #########################################################

        fig, ax, nmjd = plot_lc(r, unit='mag')
        fn = os.path.join(fermi_lightcurves_plot_dir, f"{r.name[0]}_{r['4FGL Name'].replace(' ', '_')}mag.pdf")
        print(f"saving under {fn}")
        fig.subplots_adjust(bottom=0.2)
        fig.subplots_adjust(left=0.15)
        fig.savefig(fn)
        plt.close()

        fig, ax, nmjd = plot_lc(r, unit='mag', timerange=[-100, 100])
        fn = os.path.join(fermi_lightcurves_plot_dir, f"{r.name[0]}_{r['4FGL Name'].replace(' ', '_')}mag_closeup.pdf")
        fig.subplots_adjust(bottom=0.2)
        fig.subplots_adjust(left=0.15)
        print(f"saving under {fn}")
        fig.savefig(fn)
        plt.close()

        #########################################################
        #           Plot PKS 1502+106
        #########################################################

        fgl_name = '4FGL J1504.4+1029 '
        r = fermi_matches.loc[fermi_matches['4FGL Name'] == fgl_name].iloc[0]
        fig, ax, nmjd = plot_lc(r, unit='mag', timerange=[-400, 400])
        fn = os.path.join(fermi_lightcurves_plot_dir,
                          f"{r.name[0]}_{r['4FGL Name'].replace(' ', '_')}mag_closeup.pdf")
        fig.subplots_adjust(bottom=0.2)
        fig.subplots_adjust(left=0.15)
        print(f"saving under {fn}")
        ax.set_xticks([58400, 58600, 58800, 59000])
        fig.savefig(fn)