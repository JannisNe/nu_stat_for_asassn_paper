import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.time import Time

from asassn_nu.alerts import get_alerts
from asassn_nu.style import big_fontsize, base_width, base_height, dpi, CB_color_cycle
from asassn_nu.info import output_folder_figures


def make_obs_plots():

    #########################################################
    #           Loading Alerts
    #########################################################

    alerts = get_alerts()

    Nalerts = len(alerts.Event.unique())
    N_retracted = len(alerts[alerts.retracted].Event.unique())
    Nalerts_observed = len(alerts[alerts.observed].Event.unique())
    Nalerts_unobserved = len(alerts[~alerts.observed & ~alerts.retracted & ~alerts["direction uncertainty missing"]].Event.unique())
    N_observed_2h = len(alerts[alerts["2h coverage"]>0].Event.unique())
    N_dir_missing = len(alerts[alerts["direction uncertainty missing"]])
    N_not_observable = 14
    N_observable_not_observed = 1

    print(
        f"{Nalerts} alerts in total\n"
        f"{N_retracted} alerts retracted\n"
        f"{Nalerts_observed} observed by ASAS-SN\n"
        f"{Nalerts_unobserved} not observed\n"
        f"{N_observed_2h} observed within 2 hours\n"
        f"{Nalerts_observed-N_observed_2h} observed within 14 days\n"
        f"{N_dir_missing} alerts missing direction uncertainties"
    )

    #########################################################
    #           Pie Chart Plot
    #########################################################

    labels = ['retracted', 'direction\nuncertainties\nmissing', 'observed', 'not observable', 'missed']
    sizes = [N_retracted, N_dir_missing, Nalerts_observed, N_not_observable, N_observable_not_observed]
    explode = [0, 0, 0.05, 0, 0]
    labels_with_numbers = [f"{l} ({n})" if n < 3 else l for l, n in zip(labels, sizes)]

    assert len(labels) == len(sizes)

    def make_autopct(values):
        def my_autopct(pct):
            total = sum(values)
            val = int(round(pct*total/100.0))
            return f"{val:.0f}" if val > 3 else ""
        return my_autopct

    fig, ax = plt.subplots(figsize=(base_width, base_height), dpi=dpi)
    _, _, autotxts = ax.pie(
        sizes,
        labels=labels_with_numbers,
        autopct=make_autopct(sizes),
        shadow=False,
        startangle=135,
        pctdistance=0.85,
        textprops={'fontsize': big_fontsize},
        colors=np.array(CB_color_cycle)[[6, 8, 5, 1, 7]],
        explode=explode
    )

    for autotext in autotxts:
        autotext.set_color('white')
    ax.axis('equal')
    fig.tight_layout()

    fn = os.path.join(output_folder_figures, "alert_stats_piechart.pdf")
    print("saving under", fn)
    fig.savefig(fn)

    plt.show()
    plt.close()

    #########################################################
    #           Observation Times Plot
    #########################################################

    times = Time(list(alerts["arrival time [UT]"]))
    five_st_time = Time("2019-06-01")
    alerts["five_stations"] = times > five_st_time

    two_st_yr = (five_st_time - min(times)).to("yr").value
    total_yr = (max(times) - min(times)).to("yr").value
    five_st_yr = total_yr - two_st_yr

    two_st_alerts = len(times[(times < five_st_time) & ~alerts.Event.duplicated()])
    five_st_alerts = len(times[(times > five_st_time) & ~alerts.Event.duplicated()])
    total_alerts = len(alerts[~alerts.Event.duplicated()])

    fig, (ax, ratioax) = plt.subplots(nrows=2, figsize=(base_width, 2 * base_height), dpi=dpi,
                                      gridspec_kw={'hspace': 0, 'height_ratios': [1, 1]}, sharex='col')
    fig.subplots_adjust(wspace=0.)

    unmasked = [True] * len(alerts)

    for m, c, ls, l, normto in zip(
            [unmasked, ~alerts.five_stations, alerts.five_stations],
            [CB_color_cycle[0], CB_color_cycle[1], CB_color_cycle[2]],
            ["-", ":", "--"],
            ["total", "two stations", "five stations"],
            [total_alerts, two_st_alerts, five_st_alerts]
    ):
        mbefore = (~alerts.Before.isna()) & ~alerts.Event.duplicated() & m
        sorted_before = alerts[mbefore].sort_values(by="Before", ignore_index=True, ascending=True)

        mafter = (~alerts.After.isna()) & ~alerts.Event.duplicated() & m
        sorted_after = alerts[mafter].sort_values(by="After", ignore_index=True)

        times_before = list(-sorted_before.Before[::-1]) + [0]
        times_after = [0] + list(sorted_after.After)
        N_before = list(sorted_before.index[::-1] + 1) + [0]
        N_after = [0] + list(sorted_after.index + 1)

        ax.plot(times_before, N_before, drawstyle='steps-pre', color=c, ls=ls, label=l)
        ax.plot(times_after, N_after, drawstyle='steps-post', color=c, ls=ls)

        ratioax.plot(times_before, np.array(N_before) / normto, drawstyle='steps-pre', color=c, ls=ls)
        ratioax.plot(times_after, np.array(N_after) / normto, drawstyle='steps-post', color=c, ls=ls)

    ax.set_xscale("symlog", linthresh=1)

    for a in [ax, ratioax]:
        a.axvline(0, ls='--', color='k')
        a.grid()

    minor_ticks = np.array(list(np.linspace(10, 50, 5) / 60) + list(range(1, 24)) + list(np.array(range(1, 7)) * 24) + list(
        np.array(range(1, 3)) * 24 * 7))
    xticks = [1, 24, 168]
    xticklabels = ['1h', '1d', '1w']

    ax.set_xticks(list(-np.array(xticks)) + [0] + xticks)
    ax.set_xticklabels(["-" + xtl for xtl in xticklabels] + [0] + ["" + xtl for xtl in xticklabels])
    ax.set_xticks(list(-minor_ticks) + list(minor_ticks), minor=True)
    ax.set_xticklabels([''] * len(minor_ticks) * 2, minor=True)

    ax.set_ylim([0, 60])
    ratioax.set_ylim([0, 1])
    ratioax.set_yticks(np.linspace(0, 0.8, 5))
    ax.set_xlim([min(-alerts.Before), max(alerts.After)])

    ratioax.set_xlabel('Time to neutrino alert')
    ax.set_ylabel('observed alerts')
    ratioax.set_ylabel("probability of follow-up")
    ax.legend(bbox_to_anchor=(0.32, 0.98), ncol=1, loc='upper center')

    fig.tight_layout()

    fn = os.path.join(output_folder_figures, "time_of_observation.pdf")
    print("saving under", fn)
    fig.savefig(fn)

    plt.show()
    plt.close()


if __name__ == '__main__':
    make_obs_plots()