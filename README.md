[![DOI](https://zenodo.org/badge/203128425.svg)](https://zenodo.org/badge/latestdoi/203128425) [![MIT Licence](https://img.shields.io/badge/License-MIT-red)](https://opensource.org/licenses/mit-license.php) ![alt text](https://img.shields.io/badge/Python-3.x-brightgreen 'Supported platform') ![alt text](https://img.shields.io/badge/Tested%20on-Mac%2FLinux%2FWindows-brightgreen 'Supported platform') [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jonas-ka/pi-icr-analysis/HEAD)



# PI-ICR analysis

Created on 17 July 2019 for the ISOLTRAP experiment
- V1.1 (24 June 2020): Maximum likelihood estimation was simplified based on SciPy PDF's and the CERN-ROOT6 minimizer via the iminuit package (→ great performance)

@author: Jonas Karthein<br>
@contact: karthein@mit.edu<br>
@license: MIT license

### References
[1]: https://doi.org/10.1007/s00340-013-5621-0
[2]: https://doi.org/10.1103/PhysRevLett.110.082501
[3]: https://doi.org/10.1007/s10751-019-1601-z
[4]: https://doi.org/10.1103/PhysRevLett.124.092502

[1] S. Eliseev, _et al._ Appl. Phys. B (2014) 114: 107.<br>
[2] S. Eliseev, _et al._ Phys. Rev. Lett. 110, 082501 (2013).<br>
[3] J. Karthein, _et al._ Hyperfine Interact (2019) 240: 61.<br>

### Application

The code was used to analyse data for the following publications:

[3] J. Karthein, _et al._ Hyperfine Interact (2019) 240: 61.<br>
[4] V. Manea and J. Karthein, _et al._ Phys. Rev. Lett. 124, 092502 (2020)<br>
[5] M. Mougeot, _et al._ in preparation (2020)<br>

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
