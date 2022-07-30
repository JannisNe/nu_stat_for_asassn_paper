import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerTuple
from astropy import units as u

from asassn_nu.info import output_folder_figures
from asassn_nu.style import base_height, base_width, dpi, bandcols, lc_errbar_kw, lc_uplim_kw, bandmark


#########################################################
#           3HSPJ095507+355101
#########################################################
ra, dec = "09 55 07.8818398488", "+35 51 00.885435048"
hsp_target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

# IceCube-2001072A
# Date: 20/01/07
# Time: 09:42:18.36  UT
icecube_loc = EarthLocation(lon=0, lat=-90, height=2450)
hsp_neutrino_time = Time('2020-01-07T09:42:18.36', scale='utc', location=icecube_loc)
# heliocentric correction
ltt_helio = hsp_neutrino_time.light_travel_time(hsp_target, 'heliocentric')
neutrino_time_helio = Time(hsp_neutrino_time.utc - ltt_helio, format='jd', scale='utc', location=icecube_loc)
neutrino_hjd_hsp = neutrino_time_helio.jd


# ASAS-SN data
hsp_fns = {
    "difference phot. + reference": "data/3HSPJ095507+355101_diff_plus_ref_phot.csv",
    "difference photometry": "data/3HSPJ095507+355101_diff_phot.csv",
    "aperture phot.": "data/3HSPJ095507+355101_apperture_phot.csv"
}
hsp_datas = {k: pd.read_csv(v) for k, v in hsp_fns.items()}

#########################################################
#           AT2019fdr
#########################################################
ra, dec = 257.278579, 26.855694
tywin_target = SkyCoord(ra, dec, unit=(u.deg, u.deg))

# IceCube-IC200530A
# Date: 20/01/07
# Time: 09:42:18.36  UT
tyein_neutrino_time = Time("2020-05-30T07:54:29.43", scale='utc', location=icecube_loc)
# heliocentric correction
ltt_helio = tyein_neutrino_time.light_travel_time(tywin_target, 'heliocentric')
neutrino_time_helio = Time(tyein_neutrino_time.utc - ltt_helio, format='jd', scale='utc', location=icecube_loc)
tywin_neutrino_hjd = neutrino_time_helio.jd

# ASAS-SN data
tywin_fns = {
    "difference photometry": "data/at2019fdr_asassn_dif_phot.csv",
}
tywin_datas = {k: pd.read_csv(v) for k, v in tywin_fns.items()}


#########################################################
#           Make Plot
#########################################################

def plot_missed_candidates():
    time_from_nu = 300
    k = "difference photometry"

    fig, axs = plt.subplots(ncols=2, figsize=(base_width*2, base_height), dpi=dpi, sharey='row')

    for ax, data, neutrino_time, nuname, tname in zip(
        axs,
        [hsp_datas[k], tywin_datas[k]],
        [hsp_neutrino_time, tyein_neutrino_time],
        ["IC2001072A", "IC200530A"],
        ["3HSPJ095507+355101", "AT2019fdr"]
    ):
        print(tname)
        print("closest obs before", min(neutrino_time.jd - data.HJD[data.HJD < neutrino_time.jd]))
        print("closest obs after", min(data.HJD[data.HJD > neutrino_time.jd] - neutrino_time.jd))
        yr_priot = min(data.HJD - neutrino_time.jd) / 365
        print(f"data goes to {yr_priot:.2f} years before nu")

        em = data.flux_err < 99.

        tm = abs(data["HJD"] - neutrino_time.jd) <= time_from_nu
        um = data["mag_err"] == 99.99
        print(f"{sum(um)} ")

        line = ax.axvline(neutrino_time.mjd)
        conts = []
        labels = []
        ax.set_title(tname)

        for i, b in enumerate(data.Filter.unique()):
            bm = data.Filter == b
            d_meas = data[bm & em & tm & ~um]
            d_ulim = data[bm & em & tm & um]

            a1 = ax.errorbar(
                x=Time(d_meas['HJD'], format='jd').mjd,
                y=d_meas['mag'].astype(float),
                yerr=d_meas['mag_err'],
                color=bandcols[b],
                marker=bandmark[b],
                **lc_errbar_kw
            )

            a2 = ax.plot(
                Time(d_ulim["HJD"], format="jd").mjd,
                d_ulim["Limit"],
                color=bandcols[b],
                **lc_uplim_kw
            )

            if len(d_meas):
                conts.append((a1, a2[0]))
                labels.append(f"{b}, upper limits")

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ylim = axs[0].get_ylim()
    axs[0].set_ylim(ylim[1], ylim[0])

    fig.supxlabel("MJD")
    fig.supylabel("Magnitude")
    fig.legend(
        [line] + conts, ["neutrino"] + labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.95),
        bbox_transform=fig.transFigure
    )
    fig.subplots_adjust(bottom=0.2)
    fig.subplots_adjust(top=0.7)
    fig.subplots_adjust(left=0.1)
    # fig.tight_layout()
    fig_fn = f"{output_folder_figures}at2019fdr_3hsp{k.replace(' ', '').replace('.', '')}.pdf"
    print(fig_fn)
    fig.savefig(fig_fn)

    plt.show()
    plt.close()


if __name__ == '__main__':
    plot_missed_candidates()