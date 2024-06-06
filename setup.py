import setuptools


if __name__ == '__main__':
    setuptools.setup(
        name="asassn_nu",
        version="0.1",
        author="Jannis Necker",
        author_email="jannis.necker@gmail.com",
        description="Code to produce results published  in https://arxiv.org/abs/2204.00500",
        license="MIT",
        keywords="neutrino astronomy high-energy optical",
        url="https://github.com/JannisNe/nu_stat_for_asassn_paper",
        project_urls={
            "Bug Tracker": "https://github.com/JannisNe/nu_stat_for_asassn_paper/issues",
        },
        packages=setuptools.find_packages(),
        classifiers=[
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3.10",
        ],
        python_requires='>=3.9',
        install_requires=[
            "matplotlib==3.5.2",
            "pandas==1.4.1",
            "jupyterlab==3.2.9",
            "seaborn==0.11.2",
            "numpy==1.23.1",
            "astropy==5.0.4",
            "scipy==1.9.0",
            "tqdm==4.64.0 ",
            "flarestack==2.4.5"
        ],
        package_data={'timewise': [
            'wise_flux_conversion_correction.dat'
        ]}
    )
