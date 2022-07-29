import numpy as np
import os
from astropy.coordinates import SkyCoord, get_sun
from alerts import get_alerts
from astropy.time import Time

from asassn_nu.info import output_folder_latex


def make_latex_tables():

    #########################################################
    #           Loading Alerts
    #########################################################

    alerts = get_alerts()

    #########################################################
    #           Tables
    #########################################################

    # ---------------        Observed Alerts        --------------- #

    m = alerts['observed'] & (~alerts.retracted)

    text = ""

    tot_area = 0.

    for index, row in alerts[m].iterrows():

        if (row.Event == 'IC160731A') & (row.Class == 'HESE'):
            continue

        name = str(row["Event"].lower())
        ras = row["RA Unc (rectangle) float"]
        decs = row["Dec Unc (rectangle) float"]

        if row['Event'] in ['IC190504A', 'IC200227A']:
            area = row['initial Error90 [arcmin]'] ** 2 * 4 * np.pi
        try:
            if row['Event'] == 'IC190629A':
                delta_r = 2 * np.pi

            else:
                delta_r = ras[0] - ras[1]
            delta_d = decs[0] - decs[1]

            area = delta_r * delta_d * np.cos(np.radians(float(row["Dec"])))
        except TypeError:
            print(row, '\n')
            area = np.nan

        if not np.isnan(area):
            tot_area += area * row["14d coverage"] / 100

        if np.isnan(float(row["Signalness"])):
            s = "-"
        else:
            s = f'{100. * row["Signalness"]:.0f}'

        if row.Class == 'EXTRA':
            s += "$^{*}$"

        ref = '-' if row.Event == "IC160731A" else f"\cite{{{name}}}"

        text += f'{row["Event"]} & {row["RA"]:.2f} & ${row["Dec"]:+.2f}$ & {area:.1f} & {row["24h coverage"]:.1f} & {row["14d coverage"]:.1f} & {s} & {ref} \\\ \n'.replace(
            "-", "-")

    print(text)
    fn = os.path.join(output_folder_latex, 'alert_table_observed.tex')
    print("saving under", fn)
    with open(fn, 'w') as f:
        f.write(text)

    # ---------------        Unobserved Alerts        --------------- #

    m = ~alerts['observed'] & (~alerts['retracted'])
    text = ""

    for index, row in alerts[m].iterrows():

        name = str(row["Event"].lower())
        ras = row["RA Unc (rectangle) float"]
        decs = row["Dec Unc (rectangle) float"]

        try:
            if row['Event'] in ['IC190629A', 'IC190504A', 'IC200227A']:
                continue
            else:
                delta_r = ras[0] - ras[1]
            delta_d = decs[0] - decs[1]

            area = delta_r * delta_d * np.cos(np.radians(float(row["Dec"])))
        except TypeError:
            print(row, '\n')
            area = np.nan

        sun_pos = get_sun(Time(row['arrival time [UT]']))
        alert_pos = SkyCoord(row["RA"], row["Dec"], unit='deg').transform_to('gcrs')

        text += (
            f'\t {row["Event"]} & '
            f'{row["RA"]} & '
            f'${row["Dec"]:+.2f}$ & '
            f'{area:.1f} & '
            f'{row["reason"]} & '
            f'\cite{{{name}}} \\\ \n'
        )

    print(text)

    fn = os.path.join(output_folder_latex, "alert_table_not_observed.tex")
    print("saving under", fn)
    with open(fn, 'w') as f:
        f.write(text)


if __name__ == '__main__':
    make_latex_tables()
