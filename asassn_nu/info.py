import os


output_dir = os.getenv("ASASSN_NU_OUTPUT", "../../")
output_folder_figures = os.path.join(output_dir, "figures")
output_folder_latex = os.path.join(output_dir, "latex")

data_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), "../../data"))


def make_directories():
    for d in [output_dir, output_folder_latex, output_folder_figures]:
        if not os.path.isdir(d):
            os.mkdir(d)
