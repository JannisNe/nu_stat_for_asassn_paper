import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import copy
from collections import OrderedDict
from astropy.cosmology import Planck18 as cosmo
from astropy.time import Time
import pandas as pd
from scipy.interpolate import interp1d
from flarestack.cosmo import define_cosmology_functions

from asassn_nu.alerts import get_alerts
from asassn_nu.info import output_folder_figures
from asassn_nu.style import base_width, base_height, dpi, CB_color_cycle
from asassn_nu.limits import rates, max_dl, cl, pval


# the list elements: declination range, completeness magnitude, color, linestyle
asassn_comp_mag = 17.5
asassn_color = CB_color_cycle[1]
instruments = {
    'ZTF': [(-30, 90), 20.8, CB_color_cycle[5], "-."],
    'LSST': [(-90, 10), 24, CB_color_cycle[4], "--"],
    'All-Sky LSST': [(-90, 90), 24, CB_color_cycle[6], ":"],
    'ASAS-SN': [(-90, 90), asassn_comp_mag, asassn_color, "-"]
}

ts = ['24h', '14d']


#########################################################
#           Calculate Population CDFs
#########################################################

def calculate_population_cdf():

    nsteps = 1e4
    zrange, step = np.linspace(0.0, 800, int(nsteps + 1), retstep=True)
    cdf_mpc = dict()
    zplot = 0.5 * (zrange[1:] + zrange[:-1])
    dlum_plot = cosmo.luminosity_distance(zplot)

    for label, rate in rates:
        rate_per_z, nu_flux_per_z, nu_flux_per_source, cumulative_nu_flux = \
            define_cosmology_functions(rate, 1., gamma=2.0)

        y = [0] + [x.value for x in cumulative_nu_flux(zplot)]
        y = np.array(y) / y[-1]
        dls = [0.] + [dl.value for dl in dlum_plot[1:]]
        cdf_mpc[label] = interp1d(dls, [0.] + list(y))

    return cdf_mpc


#########################################################
#           Calculate Upper Limits
#########################################################


def calculate_upper_limits(the_alerts, declination_selection=True):
    uls = dict()
    for inst, info in instruments.items():

        print(f'---------- {inst} -----------')

        if declination_selection:
            sample = the_alerts[(the_alerts.Dec > min(info[0])) & (the_alerts.Dec < max(info[0]))]
        else:
            sample = the_alerts

        print(len(sample), 'alerts used')

        median_coverages = {}

        def P_det(f, t):
            if t not in median_coverages:
                median_coverages[t] = np.median(sample[f"{t} coverage"][~np.isnan(sample[f"{t} coverage"])])
                print(t, 'median coverage', median_coverages[t], '%')

            P = 1
            for i, r in sample.iterrows():
                if r.Event == 'IC160731A' and r.Class == 'HESE':
                    continue
                P *= 1 - r.Signalness * median_coverages[t] / 100 * f
            return P

        fs = np.linspace(0, 1, 1000)
        Ps = {t: [P_det(f, t) for f in fs] for t in ts}

        inverse_pdet = {t: interp1d(iPs, fs) for t, iPs in Ps.items()}

        ul_sample = {t: ip(pval) for t, ip in inverse_pdet.items()}
        for t, iul in ul_sample.items():
            print(f'{t}: Only {iul * 100:.1f}% ({cl * 100:.0f}% CL) can lie above our limiting magnitude.')
            fig, ax = plt.subplots()
            ax.plot(fs, Ps[t], color='k')
            plt.axhline(pval)
            plt.axvline(ul_sample[t])
            ax.set_xlim((0, 1))
            ax.set_title(t)
            plt.show()
            plt.close()

        uls[inst] = ul_sample

    return uls


#########################################################
#           Make Plots
#########################################################


def make_upper_limits_plots():

    cdf_mpc = calculate_population_cdf()

    #########################################################
    #           Load Alerts
    #########################################################

    alerts = get_alerts()
    m = (~alerts.retracted & alerts.observed)
    use_alerts = alerts[m]
    hese_m = use_alerts.Class == 'HESE'
    use_alerts.loc[hese_m, 'Signalness'] = 0.

    # mimic ZTF period
    time_mask = Time(list(use_alerts["arrival time [UT]"])) > Time("2018-03-20")
    use_alert_time_masked = use_alerts[time_mask]

    # select ZTF observed events
    obs_evs = pd.read_csv("data/nu_alerts_observed.csv", skiprows=[0, 1, 2]).iloc[:-3].Event
    ztf_mask = list()
    for i, r in alerts.iterrows():
        ztf_mask.append(r.Event in list(obs_evs))

    use_alerts_ztf = alerts[ztf_mask]

    add = pd.DataFrame({
        "Event": ["IC210922A"],
        "Class": ["GOLD"],
        "RA": [60.73],
        "RA Unc (rectangle)": [[0.96, -0.66]],
        "Dec": [-4.18],
        "Dec Unc (rectangle)": [[0.42, -0.55]],
        "Signalness": [0.925]
    })

    use_alerts_ztf = pd.concat([use_alerts_ztf, add])

    N = len(use_alerts.Event.unique())
    N_time = len(use_alert_time_masked.Event.unique())
    N_ztf = len(use_alerts_ztf)
    print(N, 'alerts used')
    print(N_time, "alerts used when mimicing ZTF time")
    print(N_ztf, "ZTF alerts")

    #########################################################
    #           Calculate Upper Limits
    #########################################################

    uls = calculate_upper_limits(use_alerts)
    uls_time_masked = calculate_upper_limits(use_alert_time_masked)
    uls_ztfalerts = calculate_upper_limits(use_alerts_ztf, declination_selection=False)

    uls_assasn_res = copy.deepcopy(uls)

    with open('data/asassn_ul.pkl', 'rb') as f:
        uls_assasn_res['ASAS-SN'] = pickle.load(f)

    #########################################################
    #       Project ASAS-SN Obs on other instruments
    #########################################################

    f = cdf_mpc["SFR"]
    t = '14d'

    fig, ax = plt.subplots(
        figsize=(base_width, base_height),
        dpi=dpi,
        sharey='all',
        gridspec_kw={'wspace': 0}
    )

    abs_mag = np.linspace(-16, -27, 1000)

    for inst, d in uls.items():
        compmag = instruments[inst][1]
        color = instruments[inst][2]
        ls = instruments[inst][3]
        lw = 1

        def max_f(abs_mags):
            res = d[t] / f(max_dl(abs_mags, compmag))
            return res

        ax.plot(abs_mag, max_f(abs_mag), label=inst, color=color, ls=ls, lw=lw)
        mag_scatter = np.linspace(abs_mag[0], abs_mag[-1], 10)[1:-1]
        ax.errorbar(
            mag_scatter,
            max_f(mag_scatter),
            yerr=0.0,
            uplims=True,
            color=color,
            linestyle="",
            capsize=0,
            mew=1
        )

    ax.set_xlim((max(abs_mag), min(abs_mag) + 1))
    ax.set_ylim((0, 1))
    ax.legend()
    ax.set_xlabel('g-Band Peak Absolute Mag')
    ax.set_ylabel('$F_{L}$')
    fig.tight_layout()

    fn = os.path.join(output_folder_figures, "future_inst.pdf")
    print(fn)
    fig.savefig(fn)

    plt.show()
    plt.close()

    #########################################################
    #       Use actual ZTF results
    #########################################################

    evolution = "SFR"
    f = cdf_mpc[evolution]
    t = "14d"

    fig, ax = plt.subplots(
        figsize=(base_width, base_height),
        dpi=dpi,
        sharey='all',
        gridspec_kw={'wspace': 0}
    )

    for inst, d in uls.items():

        compmag = instruments[inst][1]
        color = instruments[inst][2]
        ls = instruments[inst][3]
        lw = 1

        if inst == 'ZTF':
            dat = pd.read_csv(f"data/ztf_nu_{evolution.lower()}.csv", sep=";", decimal=",", names=["abs_mag", "f"])
            poly_coeff = np.polyfit(dat.abs_mag, dat.f, 6)
            max_f = np.poly1d(poly_coeff)
        else:
            def max_f(abs_mags):
                res = d[t] / f(max_dl(abs_mags, compmag))
                return res

        abs_mag = np.linspace(-16, -27, 1000)
        mag_scatter = np.linspace(abs_mag[0], abs_mag[-1], 10)[1:-1]
        ax.plot(abs_mag, max_f(abs_mag), label=inst, color=color, ls=ls, lw=lw)
        ax.errorbar(
            mag_scatter,
            max_f(mag_scatter),
            yerr=0.0,
            uplims=True,
            color=color,
            linestyle="",
            capsize=0,
            mew=1
        )

    ax.set_xlim((max(abs_mag), min(abs_mag) + 1))
    ax.set_ylim((0, 1))
    ax.legend()
    ax.set_xlabel('g-Band Peak Absolute Mag')
    ax.set_ylabel('$F_{L}$')
    fig.tight_layout()
    fn = os.path.join(output_folder_figures, "future_inst_ztf_res.pdf")
    print(fn)
    fig.savefig(fn)

    plt.show()
    plt.close()

    #######################################################################
    #     Real ASAS-SN meas, real ZTF meas, LSST projected from ASAS-SN
    #######################################################################

    ordereD_uls_asassn_res = OrderedDict(uls_assasn_res)
    for k in ["ASAS-SN", "ZTF", "LSST", "All-Sky LSST"]:
        ordereD_uls_asassn_res.move_to_end(k)

    evolution = "SFR"
    f = cdf_mpc[evolution]
    t = "14d"

    fig, ax = plt.subplots(
        figsize=(base_width, base_height),
        dpi=dpi,
        sharey='all',
        gridspec_kw={'wspace': 0}
    )

    for inst, d in ordereD_uls_asassn_res.items():

        compmag = instruments[inst][1]
        color = instruments[inst][2]
        ls = instruments[inst][3]
        lw = 1

        if inst == 'ZTF':
            dat = pd.read_csv(f"data/ztf_nu_{evolution.lower()}.csv", sep=";", decimal=",", names=["abs_mag", "f"])
            poly_coeff = np.polyfit(dat.abs_mag, dat.f, 6)
            max_f = np.poly1d(poly_coeff)
        else:
            def max_f(abs_mags):
                res = d[t] / f(max_dl(abs_mags, compmag))
                return res

        abs_mag = np.linspace(-16, -27, 1000)
        mag_scatter = np.linspace(abs_mag[0], abs_mag[-1], 10)[1:-1]
        ax.plot(abs_mag, max_f(abs_mag), label=inst, color=color, ls=ls, lw=lw)
        ax.errorbar(
            mag_scatter,
            max_f(mag_scatter),
            yerr=0.0,
            uplims=True,
            color=color,
            linestyle="",
            capsize=0,
            mew=1
        )

    ax.set_xlim((max(abs_mag), min(abs_mag) + 1))
    ax.set_ylim((0, 1))
    ax.legend(
        bbox_to_anchor=(0.5, 1.01),
        loc='lower center',
        ncol=2,
        columnspacing=5.7
    )
    ax.set_xlabel('g-Band Peak Absolute Mag')
    ax.set_ylabel('$F_{L}$')
    fig.tight_layout()

    fn = os.path.join(output_folder_figures, "future_inst_ztf_res_asasn_res.pdf")
    print(fn)
    fig.savefig(fn)

    plt.show()
    plt.close()

    #######################################################################
    #     Real ZTF meas, others projected with ZTF livetime
    #######################################################################

    evolution = "SFR"
    f = cdf_mpc[evolution]
    t = "14d"

    fig, ax = plt.subplots(
        figsize=(base_width, base_height),
        dpi=dpi,
        sharey='all',
        gridspec_kw={'wspace': 0}
    )

    use_uls = copy.deepcopy(uls_ztfalerts)
    use_uls["ASAS-SN (full)"] = uls_time_masked["ASAS-SN"]

    for inst, d in use_uls.items():

        if "All-Sky" in inst:
            continue

        if inst == 'ASAS-SN (full)':
            compmag = asassn_comp_mag
            color = asassn_color
            ls = "-"
            lw = 1.5
            label = r"ASAS-SN ($N_\nu$" + f"={len(use_alert_time_masked)})"
        else:
            compmag = instruments[inst][1]
            color = instruments[inst][2]
            ls = instruments[inst][3] if inst != "ASAS-SN" else ":"
            lw = 1
            label = inst + r" ($N_\nu$" + f"={len(use_alerts_ztf)})"

        if inst == 'ZTF':
            dat = pd.read_csv(f"data/ztf_nu_{evolution.lower()}.csv", sep=";", decimal=",", names=["abs_mag", "f"])
            poly_coeff = np.polyfit(dat.abs_mag, dat.f, 6)
            max_f = np.poly1d(poly_coeff)
        else:
            def max_f(abs_mags):
                res = d[t] / f(max_dl(abs_mags, compmag))
                return res

        abs_mag = np.linspace(-16, -27, 1000)
        mag_scatter = np.linspace(abs_mag[0], abs_mag[-1], 10)[1:-1]
        ax.plot(abs_mag, max_f(abs_mag), label=label, color=color, ls=ls, lw=lw)
        ax.errorbar(
            mag_scatter,
            max_f(mag_scatter),
            yerr=0.0,
            uplims=True,
            color=color,
            linestyle="",
            capsize=0,
            mew=1
        )

    ax.set_xlim((max(abs_mag), min(abs_mag) + 1))
    ax.set_ylim((0, 1))
    ax.legend(
        bbox_to_anchor=(0.5, 1.01),
        loc='lower center',
        ncol=2,
    )
    ax.set_xlabel('g-Band Peak Absolute Mag')
    ax.set_ylabel('$F_{L}$')
    fig.tight_layout()
    fn = os.path.join(output_folder_figures, "future_inst_ztf_res_ztf_time.pdf")
    print(fn)
    fig.savefig(fn)

    plt.show()
    plt.close()


if __name__ == '__main__':
    make_upper_limits_plots()