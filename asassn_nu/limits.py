import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from scipy.interpolate import interp1d
from flarestack.cosmo import get_rate, define_cosmology_functions

from asassn_nu.alerts import get_alerts
from asassn_nu.info import output_folder_figures, data_dir
from asassn_nu.style import base_width, base_height, dpi, CB_color_cycle


upper_limits_fn = os.path.join(data_dir, "asassn_ul.pkl")
cl = 0.9
pval = 1 - cl

rates = [
    ("GRB", get_rate("GRB")),
    ("SFR", get_rate("ccsn", rate_name="strolger_15", evolution_name="strolger_15", kcc_name="strolger_15")),
]

labels = {
    "TDE": "TDE-like",
    "GRB": "GRB-like",
    "CCSN (Strolger)": "SFR-like"
}

colors = {
    "TDE": CB_color_cycle[0],
    "GRB": CB_color_cycle[1],
    "CCSN (Madau)": CB_color_cycle[3],
    "CCSN (Strolger)": CB_color_cycle[5],
    "CCSN": CB_color_cycle[5],
    "SFR": CB_color_cycle[5]
}


#########################################################
#           Calculate Upper Limits
#########################################################


def calculate_upper_limits():

    alerts = get_alerts()
    m = alerts['observed'] & (~alerts.retracted)
    observed = alerts[m]

    N = len(observed.Event.unique())
    print(N, 'alerts observed')

    hese_m = observed.Class == 'HESE'
    observed.loc[hese_m, 'Signalness'] = 0.

    def P_det(f, t):
        P = 1
        for i, r in observed.iterrows():
            if r.Event == 'IC160731A' and r.Class == 'HESE':
                continue
            P *= 1 - r.Signalness * r[f'{t} coverage'] / 100 * f
        return P

    fs = np.linspace(0, 1, 1000)
    Ps = {t: [P_det(f, t) for f in fs] for t in ['1h', '2h', '4h', '24h', '14d']}

    inverse_pdet = {t: interp1d(iPs, fs) for t, iPs in Ps.items()}

    ul = {t: ip(pval) for t, ip in inverse_pdet.items()}
    for t, iul in ul.items():
        print(f'{t}: Only {iul * 100:.1f}% ({cl * 100:.0f}% CL) can lie above our limiting magnitude.')
        fig, ax = plt.subplots()
        ax.plot(fs, Ps[t], color='k')
        plt.axhline(pval)
        plt.axvline(ul[t])
        ax.set_xlim((0, 1))
        ax.set_title(t)
        plt.show()
        plt.close()

    with open(upper_limits_fn, "wb") as f:
        pickle.dump(ul, f)


def get_upper_limits(force_new=False):
    if (not os.path.isfile(upper_limits_fn)) or force_new:
        calculate_upper_limits()

    with open(upper_limits_fn, "wb") as f:
        uls = pickle.load(f)

    return uls


def calculate_population_cdf():
    nsteps = 1e3
    zrange, step = np.linspace(0.0, 10, int(nsteps + 1), retstep=True)
    cdf_mpc = dict()
    zplot = 0.5 * (zrange[1:] + zrange[:-1])
    dlum_plot = cosmo.luminosity_distance(zplot)

    for label, rate in rates:
        rate_per_z, nu_flux_per_z, nu_flux_per_source, cumulative_nu_flux = \
            define_cosmology_functions(rate, 1., gamma=2.0)

        y = [x.value for x in cumulative_nu_flux(8.)]
        y = np.array(y) / y[-1]
        dls = [0.] + [dl.value for dl in dlum_plot]
        cdf_mpc[label] = interp1d(dls, [0.] + list(y))

    fig, ax = plt.subplots(figsize=(base_width, base_height), dpi=dpi)

    loglum_range = (1, 5)
    dlum_plot2 = np.logspace(loglum_range[0], loglum_range[1], 1000)
    for label, cdf in cdf_mpc.items():
        ax.plot(dlum_plot2, cdf(dlum_plot2), label=label, color=colors[label])

    ax.set_xscale('log')
    ax.set_ylim([0, 1])
    ax.set_xlim([10 ** loglum_range[0], 10 ** loglum_range[1]])
    ax.legend()
    ax.set_ylabel('CDF')
    ax.set_xlabel('Lumniosity Distance [Mpc]')
    fig.tight_layout()

    fn = os.path.join(output_folder_figures, "neutrino_cdfs.pdf")
    print("saving under", fn)
    fig.savefig(fn)

    plt.show()
    plt.close()

    return cdf_mpc


def max_dl(ab_mag, lim_mag):
    dl = (10. ** (0.2 * (lim_mag - ab_mag))) * (10. * u.pc)
    return dl.to(u.Mpc)


def make_upper_lim_plot():

    mag_label_map = np.array([
        ('SN IIn', 'CCSN'),
        ('TDE', 'TDE'),
        ('GRB', 'GRB'),
        ('SN IIn', 'SFR')
    ])

    lim_mag = {
        'V': {'2h': 16, '24h': 16, '14d': 16},
        'g': {'2h': 17, '24h': 17, '14d': 17}
    }

    time_to_use = {
        "GRB": ['24h', '2h'],
        "SFR": ["14d"]
    }

    labels = {
        "TDE": "TDE-like",
        "GRB": "GRB-like",
        "SFR": "SFR-like"
    }

    average_absolute_magnitudes = {
        'SN IIn': -22,
        'GRB': -27
    }

    abs_mag_band = {
        'V': np.linspace(-22, -30, 1000),
        'g': np.linspace(-22, -30, 1000)
    }

    cdf_mpc = calculate_population_cdf()
    ul = get_upper_limits()

    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(base_width * 2, base_height * 1.1), dpi=dpi, sharey='all',
                            gridspec_kw={'wspace': 0})


    for iband, (ax, band) in enumerate(zip(axs, ['V', 'g'])):
        abs_mag = abs_mag_band[band]

        for i, rate in enumerate(["SFR", "GRB"]):

            f = cdf_mpc[rate]
            for iti, ti in enumerate(time_to_use[rate]):
                iul = ul[ti]

                def max_f(abs_mags):
                    res = iul / f(max_dl(abs_mags, lim_mag[band][ti]))
                    return res

                tistr = {
                    '24h': '1 day',
                    '14d': '2 weeks',
                    '2h': '2 hours'
                }[ti]
                label = labels[rate] + ', ' + tistr
                ls = ['-', '-.'][iti]
                ax.plot(abs_mag, max_f(abs_mag), label=label, color=colors[rate], ls=ls)
                mag_scatter = np.linspace(abs_mag[0], abs_mag[-1], 10)[1:-1]
                ax.errorbar(
                    mag_scatter,
                    max_f(mag_scatter),
                    yerr=0.0,
                    uplims=True,
                    color=colors[rate],
                    linestyle="",
                    capsize=0,
                    mew=1
                )

                Mag = average_absolute_magnitudes[mag_label_map[:, 0][mag_label_map[:, 1] == rate][0]]
                print(f"{rate} {band}: {max_f(Mag) * 100:.2f}% over {Mag} at {ti}")

        Nticks = max(abs_mag) - min(abs_mag) + 1

        _s = max(abs_mag)
        _e = min(abs_mag)
        if iband == 0:
            _e += 1
            Nticks -= 1

        xticks = np.linspace(_s, _e, num=int(Nticks))
        ax.set_xticks(xticks)
        ax.set_xlim((max(abs_mag), min(abs_mag)))
        ax.set_ylim((0, 1))

        ax.set_xlabel(f'{band} band')
    axs[0].set_ylabel('$F_{L}$')
    axs[0].legend(bbox_to_anchor=(1, 1.25), ncol=3, loc='upper center')

    fig.supxlabel('Peak Absolute Mag', y=0.1)
    fig.tight_layout()

    fn = os.path.join(output_folder_figures, f"limits_{cl*100:.0f}cl_notde_separate_bands.pdf")
    print("saving under", fn)
    fig.savefig(fn)

    plt.show()
    plt.close()


if __name__ == '__main__':
    make_upper_lim_plot()