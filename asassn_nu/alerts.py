import pandas as pd
import numpy as np
import os
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u

from asassn_nu.info import data_dir


alerts_fn = os.path.join(data_dir, 'ASASSN_sample_paper_IceCube_info.csv')

notice_summary_urls = {
    'BRONZE_GOLD': "https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html",
    'EHE': "https://gcn.gsfc.nasa.gov/amon_ehe_events.html",
    'HESE': "https://gcn.gsfc.nasa.gov/amon_hese_events.html"
}


def get_summary_table(k, renew=False):
    fn = os.path.join(data_dir, f"gcn_notice_summary_{k}.html")
    if not os.path.isfile(fn) or renew:
        _t = pd.read_html(requests.get(notice_summary_urls[k]).text)[0]
        _t.to_html(fn)
    else:
        _t = pd.read_html(fn, index_col=0)[0]
    return _t


offline_alert_times = {
    "IC200107A": "2020-01-07 09:42:18.36",
    "IC210503A": "2021-05-03 22:19:32.96",
    "IC210717A": "2021-07-17 15:45:19.48"
}


def get_notice_info(ic_name, alert_class, verbose=True):
    
    if ('BRONZE' in alert_class) or ('GOLD' in alert_class):
        _alert_class = 'BRONZE_GOLD'
    else:
        _alert_class = alert_class
        
    _pos_ext = ' [deg]' if _alert_class == 'BRONZE_GOLD' else ''
    _error_ext = '90 [arcmin]' if _alert_class == 'BRONZE_GOLD' else ''
        
    summary_table = get_summary_table(_alert_class)
    
    _date_from_name = ic_name.split('IC')[-1][:-1]
    _dates_in_table = summary_table['EVENT', 'Date'].apply(lambda x: x.replace('/', ''))
    _mask = _dates_in_table == _date_from_name
    
    if 'Rev' in summary_table['EVENT'].columns:
        _mask = _mask & (summary_table['EVENT', 'Rev'] == 0)
    
    _selected = summary_table[_mask]
    
    if len(_selected) != 1:
        
        if 'IC160427A' in ic_name:
            if verbose:
                print(f'{ic_name}: selecting the third notice of {len(_selected)}')
            _ind = 1
            
        elif len(_selected) == 2:
            
            if verbose:
                print(f"{ic_name}: Two matching dates.")
                
            _ras = _selected["OBSERVATION"][f"RA{_pos_ext}"]
            _decs = _selected["OBSERVATION"][f"Dec{_pos_ext}"]
            _coords = SkyCoord(_ras, _decs, unit='deg')
            _sep = _coords[0].separation(_coords[1]).deg
            if verbose:
                print(f"\t{_sep:.2f} degrees apart")
                
            if _sep > 1:
                if verbose:
                    print(f"\tassuming it's two alerts at the same day")
                tstrings = np.array([f"20{_s['EVENT', 'Date'].replace('/','-')}T{_s['EVENT', 'Time UT']}"
                           for _, _s in _selected.iterrows()])
                times = Time(tstrings)
                _ind = np.argmin(times) if ic_name.endswith('A') else np.argmax(times) if ic_name.endswith('B') else None
                
            else: 
                if verbose:
                    print(f"\tassuming second notice is refined info from circular. choosing first one")
                _ind = 0
                
        else:
            raise Exception(f"More than one entry for {ic_name}: {_selected}")
            
    else:
        _ind = 0
        
    _selected = _selected.iloc[_ind]
    return _selected, _pos_ext, _error_ext


def parse_notice_info(ic_name, alert_class, verbose=True):
    
    if ic_name in ['IC210503A', 'IC200107A', 'IC210717A']:
        if verbose: print(f"{ic_name}: No notice becasue selected offline")
        return offline_alert_times[ic_name], None, None, None
    
    _selected, _pos_ext, _error_ext = get_notice_info(ic_name, alert_class, verbose=True)
    _date = _selected["EVENT"]["Date"].replace("/", "-")
    _obstime = _selected["EVENT"]["Time UT"]
    _ra = _selected["OBSERVATION"][f"RA{_pos_ext}"]
    _dec = _selected["OBSERVATION"][f"Dec{_pos_ext}"]
    _error90 = _selected["OBSERVATION"][f"Error{_error_ext}"]
    
    _arrivaltime = f"20{_date} {_obstime}"
    
    return _arrivaltime, _ra, _dec, _error90
    
    
hese_sig_table = pd.DataFrame([
    (6000, 1.09, 3.73),
    (6500, 1.00, 2.81),
    (7000, 0.91, 1.16),
    (7500, 0.84, 0.92),
], columns=['Charge', 'ns', 'nb'])


def hese_signalness(charge, verbose=True):
    m = hese_sig_table.Charge <= charge
    r = hese_sig_table.iloc[np.argmax(hese_sig_table[m].Charge)]
    _signalness = r.ns / (r.nb + r.ns)
    if verbose: print(f"Signalness={_signalness:.2f} for charge {charge}")
    return _signalness


def get_closest_obs_df():
    with open(os.path.join(data_dir, "image_time"), "r") as f:
        lines = f.read().split("\n")

    dat = list()

    for l in lines[1:-1]:
        name = l[:9]
        times = l[11:].split("\t")
        times_formatted = list()
        for t in times:
            if t:
                if t != "NaN":
                    val, unit = t.split(" ")
                    unit = "s" if unit == "sec" else unit
                    times_formatted.append((float(val) * u.Unit(unit)).to("h").value)
                else:
                    times_formatted.append(np.nan)
            
        dat.append(tuple([name] + times_formatted))
        
    closest_obs = pd.DataFrame(dat, columns=lines[0].replace(" ", "").split("\t")[:-1])
    return closest_obs
    

def make_alerts():
    obs = pd.read_csv(os.path.join(data_dir, "nu_alerts_observed.csv"), skiprows=[0, 1, 2])
    obs = obs[~np.isnan(obs["RA"])]
    non = pd.read_csv(os.path.join(data_dir, "nu_alerts_unobserved.csv"), skiprows=[0, 1], usecols=range(11))
    comb = pd.concat([obs, non], ignore_index=True)

    # Splitting the EHE and HESE info into two rows
    m = comb['Event'] == 'IC160731A'
    comb.loc[m, 'Class'] = 'EHE'
    to_append = comb.loc[m].copy()
    to_append['Class'] = 'HESE'
    comb = comb.append(to_append)

    new_cols = ['arrival time [UT]','initial RA','initial Dec','initial Error90 [arcmin]']
    for c in new_cols:
        comb[c] = np.nan

    for j, (i, row) in enumerate(comb.iterrows()):
        m = (comb.Event == row.Event) & (comb.Class == row.Class)
        comb.loc[m, new_cols] = parse_notice_info(row['Event'], row['Class'])
        
    comb.loc[comb.Class == 'EXTRA', 'Signalness'] = 0.5
    comb.loc[comb.Class == 'HESE', 'Signalness'] = np.nan

    comb['retracted'] = (comb['Rejection reason'] == 'Alert retraction') | (comb['Rejection reason'] == 'Alert Retraction')
    comb['direction uncertainty missing'] = False
    comb.loc[(comb.Event == 'IC190504A') | (comb.Event == 'IC200227A'), 'direction uncertainty missing'] = True
    
    comb['reason'] = ''
    for k in ['retracted', 'direction uncertainty missing']:
        comb.loc[comb[k], 'reason'] = k

    keep_cols = [
        'Event', 
        'Class', 
        'RA', 
        'RA Unc (rectangle)', 
        'Dec',
        'Dec Unc (rectangle)', 
        'arrival time [UT]', 
        'Signalness',
        'initial RA', 
        'initial Dec', 
        'initial Error90 [arcmin]',
        'retracted',
        'direction uncertainty missing',
        'reason'
    ]
    
    out = comb[keep_cols].sort_values('Event')
    print(f"saved under {alerts_fn}")
    out.to_csv(alerts_fn)
    

def get_alerts():
    if not os.path.isfile(alerts_fn):
        make_alerts()
        
    alerts = pd.read_csv(alerts_fn, index_col=0)

    def floatify(x):
        if isinstance(x, str):
            return [float(s) for s in x.replace('[', '').replace(']', '').split(',')]

    for col in ["RA Unc (rectangle)", "Dec Unc (rectangle)"]:
        alerts[col + ' float'] = alerts[col].apply(floatify)

    names = ['Event', '1h', '2h', '4h', '8h', '12h', '24h', '2d', '3d', '4d', '5d', '6d', '7d', '14d']
    for i in range(1, len(names)):
        names[i] += ' coverage'

    covered = pd.read_csv(os.path.join(data_dir, "coverage_time_Jannis"), sep='\t', names=names, skiprows=1)

    for k in names[1:]:
        alerts[k] = np.nan
        for i, r in covered.iterrows():
            m = alerts.Event == r.Event
            alerts.loc[m, k] = r[k]

    alerts['observed'] = alerts['14d coverage'] > 0
    alerts.loc[~alerts.observed & ~alerts.retracted & ~alerts["direction uncertainty missing"], 'reason'] = 'proximity to sun'
    alerts.loc[alerts.Event == 'IC200421A', 'reason'] = 'operation'
    
    # -------------------------------- ADD CLOSEST OBSERVATION -------------------------------- #
    closest_obs = get_closest_obs_df()
    
    for i, r in closest_obs.iterrows():
        alerts_m = alerts.Event == r.Event
        alerts.loc[alerts_m, "controll_Event"] = r.Event
        
        for k in ["Before", "After"]:
            alerts.loc[alerts_m, k] = r[k]

    controll_m = alerts.Event[alerts.observed] == alerts.controll_Event[alerts.observed]
    if not np.all(controll_m):
        raise Exception(f"alerts missmatch!\n{alerts[alerts.observed][controll_m].to_string()}")
        
    alerts["closest_obs"] = np.fmin(alerts.After, alerts.Before)
    alerts["closest_obs_side"] = "Before"
    alerts.loc[alerts.After <= alerts.Before, "closest_obs_side"] = "After"
    alerts.loc[alerts.closest_obs.isna(), "closest_obs_side"] = np.nan
    
    return alerts


if __name__ == "__main__":
    make_alerts()
