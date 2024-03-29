[![DOI](https://zenodo.org/badge/203128425.svg)](https://zenodo.org/badge/latestdoi/203128425) [![MIT Licence](https://img.shields.io/badge/License-MIT-red)](https://opensource.org/licenses/mit-license.php) ![alt text](https://img.shields.io/badge/Python-3.x-brightgreen 'Supported platform') ![alt text](https://img.shields.io/badge/Tested%20on-Mac%2FLinux%2FWindows-brightgreen 'Supported platform') [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jonas-ka/pi-icr-analysis/HEAD)



# PI-ICR analysis

Created on 17 July 2019 for the ISOLTRAP experiment
- V1.1 (24 June 2020): Maximum likelihood estimation was simplified based on SciPy PDF's and the CERN-ROOT6 minimizer via the iminuit package (→ great performance)
- V1.2 (20 February 2021): Preparations for scientific publication and iminuit v2 update integration

@author: Jonas Karthein<br>
@contact: karthein@mit.edu<br>
@license: MIT license

### Scientific publication on this work
[1]: https://doi.org/10.1016/j.cpc.2021.108070

Please read the following article for details on the calculation and usage:<br>
[1] J. Karthein, _et al._ Comp. Phys. Comm. 267, 108070 (2021).<br>

### References
[2]: https://doi.org/10.1007/s00340-013-5621-0
[3]: https://doi.org/10.1103/PhysRevLett.110.082501
[4]: https://doi.org/10.1007/s10751-019-1601-z
[5]: https://doi.org/10.1103/PhysRevLett.124.092502
[6]: https://doi.org/10.1038/s41567-021-01326-9

[2] S. Eliseev, _et al._ Appl. Phys. B (2014) 114: 107.<br>
[3] S. Eliseev, _et al._ Phys. Rev. Lett. 110, 082501 (2013).<br>
[4] J. Karthein, _et al._ Hyperfine Interact (2019) 240: 61.<br>

### Application

The code was used to analyse data for the following publications:

[4] J. Karthein, _et al._ Hyperfine Interact (2019) 240: 61.<br>
[5] V. Manea & J. Karthein, _et al._ Phys. Rev. Lett. 124, 092502 (2020)<br>
[6] M. Mougeot, D. Atanasov, J. Karthein, _et al._ Nature Physics 17, 1099 (2021)<br>

### Introduction

The following code was written to reconstruct raw Phase-Imaging Ion-Cyclotron-Resonance (PI-ICR) data, to fit PI-ICR position information and calculate a frequency using the patter 1/2 scheme described in Ref. [1] and to determine a frequency ratio between a measurement ion and a reference ion. Additionally, the code allows to analyze isomeric states separated in pattern 2.
 data, to fit PI-ICR position information and calculate a frequency using the patter 1/2 scheme described in Ref. [1] and to determine a frequency ratio between a measurement ion and a reference ion. Additionally, the code allows to analyze isomeric states separated in pattern 2.

### Dependencies

* [pandas](https://pandas.pydata.org/) (data storage and calculation)
* [NumPy](http://www.numpy.org/) (calculation)
* [SciPy](http://www.scipy.org/) (PDFs, least squares estimation)
* [Matplotlib](https://matplotlib.org/) (plotting)
* [configparser](https://docs.python.org/3/library/configparser.html) (configuration file processing)
* [Jupyter](https://jupyter.org) (Python notebook environment)
* [iminuit](https://iminuit.readthedocs.io/en/stable/) (CERN-ROOT6 minimizer)

All packages can be fetched using pip:

`pip install pandas numpy matplotlib scipy configparser jupyter iminuit`

### Installation and execution

The whole analysis code is based on the main jupyter notebook. It can be run locally using Python 3.7+ by installing all dependencies using the `pip` command from the section above. Next, navigate to this package's folder and run `jupyter notebook pi-icr-analysis.ipynb`. Alternatively, the code can be run in the cloud using services like CERN's SWAN platform, Google Colab or Binder (use button on top of this README). The analysis code provides a full set of example data to run the entire analysis and is already set up to do so. If the execution fails in the "Frequency-ratio calculation" section, please modify the line 129 in `bin/freq_ratio.py` to one of the following values: [-10, -1, -0.001, 0.001, 1, 10]. This varys the starting parameters of the frequency ratio fit, which might get stuck for such a small set of (example) data. In case of a real experiment with many more data points, this does not pose any issues. Have fun :-)

### Folder structure

The analysis code depends on two local packages for data reconstruction and frequency-ratio fitting. They can both be found in the `bin` folder. There, ISOLTRAP's Windows conversion executable for the binary raw data from our DLD40 detector to plain text files is also located (not required for external experiments). The example data with their full analysis and plots can be found in the `data` folder and subfolders for the ion of interest (`ioi`) and the reference (`ref`) data. The alternative fitting approaches discussed in the scientific publication including the Monte-Carlo simulations and their results/plots can be found in the `alternative-approaches` folder using the jupyter notebook `PI-ICR-phase-display-simulation.ipynb`.
