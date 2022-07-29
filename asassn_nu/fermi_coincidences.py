import matplotlib
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, requests, pickle
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy import units as u
from style import output_folder_figures, big_fontsize, base_width, base_height, dpi
import seaborn as sns
import json
from alerts import get_alerts
from astropy.time import Time


fermi_matches_fn = "data"


def make_fermi_mathces():

    #########################################################
    #           Loading Alerts
    #########################################################

    alerts = get_alerts()
    m = ~alerts.retracted

    #########################################################
    #           Loading Fermi catalogues
    #########################################################

    cols = ['3FHL Name', '4FGL Name', 'Fermi ra', 'Fermi dec', 'Counterpart', 'Counterpart ra', 'Counterpart dec',
            'Counterpart prob']
    fermi_matches = pd.DataFrame(columns=cols)

    for cat, fn in zip(['3FHL', '4FGL'], ["data/gll_psch_v13.fit", "data/gll_psc_v27.fit"]):

        with fits.open(fn) as fermi_cat:
            tab = fermi_cat[1]

        fermi_ra = tab.data["RAJ2000"]
        fermi_dec = tab.data["DEJ2000"]
        fermi_coords = SkyCoord(fermi_ra, fermi_dec, unit='deg')

        for j, (i, r) in enumerate(alerts[m].iterrows()):

            #########################################################
            #           Find matches to neutrino alerts
            #########################################################

            c = SkyCoord(r['RA'], r['Dec'], unit='deg')
            delta_ra = fermi_coords.ra - c.ra
            delta_dec = fermi_coords.dec - c.dec

            try:
                if r['Dec Unc (rectangle) float']:
                    dec_unc = Angle(r['Dec Unc (rectangle) float'] * u.deg)

                    if r.Event == "IC190629A":
                        ra_unc = Angle([360, 360] * u.deg)
                    else:
                        ra_unc = Angle(r['RA Unc (rectangle) float'] * u.deg)

                    match_mask = (delta_dec <= max(dec_unc)) & (delta_dec >= min(dec_unc)) & (delta_ra <= max(ra_unc)) & (
                                delta_ra >= min(ra_unc))

                else:
                    print(r)
                    _sep = Angle(r["initial Error90 [arcmin]"] * u.arcmin)
                    sep = c.separation(fermi_coords)
                    match_mask = sep <= _sep

            except TypeError as e:
                print(r, "\n", e)

            _matches = tab.data[match_mask]

            if cat == '4FGL':
                matches = np.array([
                    _matches["ASSOC_FHL"],
                    _matches['Source_Name'],
                    _matches['RAJ2000'],
                    _matches['DEJ2000'],
                    _matches['ASSOC1'],
                    _matches['RA_Counterpart'],
                    _matches['DEC_Counterpart'],
                    np.minimum(np.array(_matches['ASSOC_PROB_BAY']).astype(float),
                               np.array(_matches['ASSOC_PROB_LR']).astype(float))
                ]).T
                columns = cols

            else:
                matches = np.array([
                    _matches['Source_Name'],
                    _matches['RAJ2000'],
                    _matches['DEJ2000'],
                    _matches['ASSOC1']
                ]).T
                columns = [f"{cat} Name", "Fermi ra", "Fermi dec", 'Counterpart']

            if len(_matches) > 0:
                print(f"{len(_matches)} matches in {cat} for {r.Event}")
                l = np.array([(r.Event, cat, j) for j in range(len(_matches))]).T
                index = pd.MultiIndex.from_arrays(l)
                ifermi_matches = pd.DataFrame(matches, columns=columns, index=index)
                fermi_matches = fermi_matches.append(ifermi_matches)

    fermi_matches = pd.DataFrame(
        fermi_matches,
        index=pd.MultiIndex.from_tuples(
            fermi_matches.index,
            names=['Event', 'Fermi Cat', 'Macth ind']
        )
    )

    #######################################################################################
    #   Match the 3FHL and the 4FGL Catalogues and remove the duplicat entries
    #######################################################################################

    drop_m = fermi_matches['3FHL Name'].duplicated(keep=False) & fermi_matches['4FGL Name'].isna()
    print(f"found {len(fermi_matches[drop_m])} duplicates")
    fermi_matches = fermi_matches[~drop_m]
    print(len(fermi_matches))
    print(len(fermi_matches[fermi_matches["4FGL Name"].isna() & (~fermi_matches["3FHL Name"].isna())]))

    for x in ['ra', 'dec']:
        counterpart = fermi_matches[f"Counterpart {x}"]
        fermi = fermi_matches[f"Fermi {x}"]
        xm = (~counterpart.isna()) & (np.array(fermi_matches['Counterpart prob']).astype(float) > 0)
        fermi_matches[x] = fermi
        fermi_matches.loc[xm, x] = counterpart[xm]

    has_cp = ~np.isnan(np.array(fermi_matches['Counterpart dec']).astype(float))

    for k in ['Counterpart', 'Counterpart ra', 'Counterpart dec']:
        fermi_matches.loc[~has_cp, k] = ""

    ncp = len(fermi_matches[has_cp])
    print(f"{ncp} with counterpart, {len(fermi_matches) - ncp} without")
    fermi_matches.to_csv('data/fermi_matches.csv')


def make_fermi_matches_plots():

    #########################################################
    #           Loading Alerts
    #########################################################

    alerts = get_alerts()
    m = ~alerts.retracted

    cmap = plt.get_cmap()

    for i, r in alerts.iterrows():

        try:
            sources = fermi_matches.loc[pd.IndexSlice[r.Event, :, :]]
        except KeyError:
            continue

        sources = sources[~sources['4FGL Name'].isna()]

        if len(sources) > 0:
            print(sources.to_string(index=False))

            fig, ax = plt.subplots(figsize=(base_width, base_height), dpi=dpi)

            ax.set_aspect(1, adjustable='datalim')
            ax.scatter(r['RA'], r['Dec'], color=cmap(0), marker='s', s=2)
            if r['RA Unc (rectangle) float'] and r['Dec Unc (rectangle) float']:
                left_bottom_corner = (
                r['RA'] + r['RA Unc (rectangle) float'][1], r['Dec'] + r['Dec Unc (rectangle) float'][1])
                extent = (r['RA Unc (rectangle) float'][0] - r['RA Unc (rectangle) float'][1],
                          r['Dec Unc (rectangle) float'][0] - r['Dec Unc (rectangle) float'][1])
                final_e = plt.Rectangle(left_bottom_corner, extent[0], extent[1], fc=cmap(0), alpha=0.4, ec=None, ls='')
                ax.add_patch(final_e)
            else:
                print(r)

            ax.scatter(r['initial RA'], r['initial Dec'], color=cmap(0.2), s=2)
            initial_e = r['initial Error90 [arcmin]'] * u.arcmin
            initial_contour = plt.Circle((r['initial RA'], r['initial Dec']), initial_e.to('deg').value, alpha=0.4,
                                         fc=cmap(0.2), ec=None, ls='')
            ax.add_patch(initial_contour)

            ax.scatter(np.array(sources['ra']).astype(float),
                       np.array(sources['dec']).astype(float), marker='x', color=cmap(0.8))
            fgl = ~sources['4FGL Name'].isna()
            for _, s in sources[fgl].iterrows():
                ylim = ax.get_ylim()
                yext = abs(ylim[1] - ylim[0])
                n = s['Counterpart'] if s['Counterpart'] else s['4FGL Name']
                ax.annotate(n, (float(s.ra), float(s.dec) + yext / 12), ha='center', color=cmap(0.8), fontsize=6,
                            snap=True)

            ax.set_title(r.Event)
            major_base = 20 if r.Event == "IC190331A" in cc else 1
            ax.xaxis.set_major_locator(MultipleLocator(base=major_base))
            ax.yaxis.set_major_locator(MultipleLocator(base=major_base))
            ax.grid()

            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')

            fig.tight_layout()
            fig.savefig(os.path.join(output_folder, f"{r.Event}.pdf"))
            try:
                plt.show()
            finally:
                plt.close()



if __name__ == '__main__':
    make_fermi_mathces()