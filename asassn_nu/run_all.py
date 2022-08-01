from asassn_nu.info import make_directories
from asassn_nu.observation_stats import make_obs_plots
from fermi_coincidences import make_fermi_matches, make_fermi_matches_plots
from asassn_nu.fermi_matches_lightcurves import plot_fermi_matches_lightcurves
from asassn_nu.candidates import make_candidates_plot, make_at2020rng_plot
from asassn_nu.limits import make_upper_lim_plot
from asassn_nu.limits_other_instruments import make_upper_limits_plots_other
from asassn_nu.missed_candidates import plot_missed_candidates
from asassn_nu.observations_by_other_telescopes import check_obs_by_others


if __name__ == '__main__':

    make_directories()
    make_obs_plots()
    plot_fermi_matches_lightcurves()
    make_candidates_plot()
    make_at2020rng_plot()
    make_upper_lim_plot()
    make_upper_limits_plots_other()
    plot_missed_candidates()
    check_obs_by_others()
