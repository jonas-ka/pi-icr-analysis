# PI-ICR analysis

Created on 17 July 2019 for the ISOLTRAP experiment

@author: Jonas Karthein<br>
@contact: jonas.karthein@cern.ch<br>
@license: MIT license

### References
[1]: https://doi.org/10.1007/s00340-013-5621-0
[2]: https://doi.org/10.1103/PhysRevLett.110.082501
[3]: https://doi.org/10.1007/s10751-019-1601-z

[1] S. Eliseev, _et al._ Appl. Phys. B (2014) 114: 107.<br>
[2] S. Eliseev, _et al._ Phys. Rev. Lett. 110, 082501 (2013).<br>
[3] J. Karthein, _et al._ Hyperfine Interact. (2019) 240: 61.<br>

### Application

The code was used to analyse data for the following publications:

[3] J. Karthein, _et al._ Hyperfine Interact (2019) 240: 61.<br>
[5] V. Manea and J. Karthein, _et al._ submitted to Phys. Rev. Lett. (2019)<br>
[6] M. Mougeot, _et al._ in preparation (2019)<br>

### Introduction

The following code was written to reconstruct raw Phase-Imaging Ion-Cyclotron-Resonance (PI-ICR) data, to fit PI-ICR position information and calculate a frequency using the patter 1/2 scheme described in Ref. [1] and to determine a frequency ratio between a measurement ion and a reference ion. Additionally, the code allows to analyze isomeric states separated in pattern 2.

### Required software and libraries

The following code was written in Python 3.7. The required libraries are listed below with a rough description for their task in the code. It doesn't claim to be a full description of the library.
* pandas (data storage and calculation)
* numpy (calculation)
* matplotlib (plotting)
* mle (maximum likelihood estimation)
* scipy (chi square fitting)
* configparser (configuration file processing)
* jupyter (Python notebook environment)

The Python-mle package can be derived from https://github.com/ibab/python-mle, all other packages can be fetched using pip:

`pip3 install pandas numpy matplotlib scipy configparser jupyter`
