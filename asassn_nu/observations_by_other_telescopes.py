import pandas as pd
import numpy as np
import os
import json
import requests
from asassn_nu.alerts import get_alerts
from asassn_nu.info import data_dir


optical_telescopes = [
    "MASTER",
    "Pan-STARRS",
    "Lick",
    "KAIT",
    "Kanata",
    "GOTO",
    "Zwicky Transient Facility",
    "AbAO",
    "DECam",
    "DDOTI"
]

gcn_circulars_fn = os.path.join(data_dir, "gcn_reports_on_observed_neutrinos.csv")
gcn_circulars_txt_fn = os.path.join(data_dir, "gcn_reports_on_observed_neutrinos.txt")


def get_alert_infos():

    if not os.path.isfile(gcn_circulars_fn):

        alerts = get_alerts()
        gcn_circulars = pd.DataFrame(columns=["subjects"])

        for n in alerts.Event.unique():
            url = f"https://gcn.gsfc.nasa.gov/other/icecube_{n[2:]}.gcn3"

            response = requests.get(url)

            if response.status_code != 200:
                print(f"{response.status_code} {response.text} for {n}")

            subjects = [l for l in response.text.split("\n") if l.startswith("SUBJECT")]
            gcn_circulars.loc[n] = [subjects]

        gcn_circulars.to_csv(gcn_circulars_fn)

    return pd.read_csv(gcn_circulars_fn, index_col=0)


def check_obs_by_others():

    gcn_circulars = get_alert_infos()
    alerts = get_alerts()

    txt = ""
    for n, r in gcn_circulars.iterrows():
        txt += f"----- {n} ------\n{json.dumps(r.subjects, indent=4)}\n\n\n"

    with open(gcn_circulars_txt_fn, "w") as f:
        f.write(txt)

    alerts["observed_by_others"] = False
    alerts["observed_by"] = None

    for n, r in gcn_circulars.iterrows():
        subjects_list = json.loads(r.subjects.replace("'", "XXX").replace('"', "'").replace("XXX", '"'))
        observed_by_mask = np.array([np.any([t in s for s in subjects_list]) for t in optical_telescopes])
        other_obs = np.any(observed_by_mask)
        alerts.loc[alerts.Event == n, "observed_by_others"] = other_obs
        alerts.loc[alerts.Event == n, "observed_by"] = json.dumps(list(np.array(optical_telescopes)[observed_by_mask]))

    only_obserevd_by_asassn_m = alerts.observed & (~alerts.observed_by_others) & (~alerts.retracted)
    only_observed_by_asassn_N = len(alerts[only_obserevd_by_asassn_m].Event.unique())
    unretracted_alerts_N = len(alerts[~alerts.retracted].Event.unique()) - 2  # two have no localisation
    fraction = only_observed_by_asassn_N / unretracted_alerts_N
    print(
        f"{only_observed_by_asassn_N} out of "
        f"{unretracted_alerts_N} alerts only obserevd by ASASSN: "
        f"{fraction * 100:.2f}%"
    )


if __name__ == '__main__':
    check_obs_by_others()
