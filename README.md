# ASAS-SN follow-up of IceCube high-energy neutrino alerts #

This repository includes code to produce the analyses published in 
[Necker et al. (2022)](https://ui.adsabs.harvard.edu/link_gateway/2022MNRAS.tmp.2251N/doi:10.1093/mnras/stac2261).

Set up your environment tun 
```bash
pip install ./
```

Specify your output directory
```bash
export ASASSN_NU_OUTPUT=</path/to/directory>
```

Produce all plots in the paper run
```bash
python asassn_nu/run_all.py
```
