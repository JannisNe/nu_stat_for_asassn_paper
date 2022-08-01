from asassn_nu.info import make_directories
from asassn_nu.observation_stats import make_obs_plots
from fermi_coincidences import make_fermi_matches, make_fermi_matches_plots
from asassn_nu.fermi_matches_lightcurves import plot_fermi_matches_lightcurves
from candidates import make_candidates_plot, make_at2020rng_plot
from limits import make_upper_lim_plot
from limits_other_instruments import make_upper_limits_plots_other
from missed_candidates import plot_missed_candidates

import shutil


if __name__ == '__main__':

    fcts = [
        make_directories,
        make_obs_plots,
        make_fermi_matches,
        make_fermi_matches_plots,
        plot_fermi_matches_lightcurves,
        make_candidates_plot,
        make_at2020rng_plot,
        make_upper_lim_plot,
        make_upper_limits_plots_other,
        plot_missed_candidates
    ]

    for f in fcts:
        while True:
            try:
                f()
                break
            except FileNotFoundError as e:
                missing_fn = str(e).replace("'", "").split(":")[-1].replace(" ", "")
                actual_fn = missing_fn.replace("data", "data_old")
                print("copying", actual_fn, "to", missing_fn)
                shutil.copy(actual_fn, missing_fn)

    # make_directories()
    # make_obs_plots()
    # plot_fermi_matches_lightcurves()
    # make_candidates_plot()
    # make_at2020rng_plot()
    # make_upper_lim_plot()
    # make_upper_limits_plots_other()
    # plot_missed_candidates()
