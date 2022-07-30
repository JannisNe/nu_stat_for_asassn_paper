# ASAS-SN follow-up of IceCube high-energy neutrino alerts #

This repository includes code to produce the analyses published in 
[Necker et al. (2022)](https://arxiv.org/abs/2204.00500).

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