import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks")
#plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{romanbar}')
plt.rcParams["font.family"] = "sans-serif"

output_folder = "../figures/"

dpi = 300

fontsize = 7.
big_fontsize = 10.0
small_fontsize = 5

golden_ratio = 1.618

base_width = 4.0
base_height = base_width/golden_ratio

margin_width = 0.5 * base_width
margin_height = margin_width/golden_ratio

full_width = 1.5 * base_width
full_height_landscape = full_width/golden_ratio
full_height_a4 = 11.75/8.25 * full_width

cmap = "rocket"

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

bandcols = {'g': CB_color_cycle[1], 'V': CB_color_cycle[0], 'r': CB_color_cycle[2]}
bandmark = {'g': 's', 'V': 'o', 'r': 'D'}

lc_errbar_kw = {'ls': '', 'markeredgecolor': 'k', 'ecolor': 'k', 'markersize': 3, 'capsize': 2, 'mew': 0.5, 'elinewidth': 0.5}
lc_uplim_kw = {'marker': 'v', 'ls': '', 'markersize': 2, 'mfc': None, 'mew': 0.5, 'fillstyle': 'none', 'zorder': 0, 'alpha': 0.5}