from asassn_nu.info import make_directories
from asassn_nu.observation_stats import make_obs_plots
from asassn_nu.fermi_matches_lightcurves import plot_fermi_matches_lightcurves
from candidates import make_candidates_plot, make_at2020rng_plot
from limits import make_upper_lim_plot
from limits_other_instruments import make_upper_limits_plots_other
from missed_candidates import plot_missed_candidates


if __name__ == '__main__':
    make_directories()
    make_obs_plots()
    plot_fermi_matches_lightcurves()
    make_candidates_plot()
    make_at2020rng_plot()
    make_upper_lim_plot()
    make_upper_limits_plots_other()
    plot_missed_candidates()