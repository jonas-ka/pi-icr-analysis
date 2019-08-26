# PI-ICR analysis
Created on 17 July 2019 for the ISOLTRAP experiment

@author: Jonas Karthein
@contact: jonas.karthein@cern.ch
@license: MIT license

### References
[1] S. Eliseev, et al. Appl. Phys. B (2014) 114: 107.
[2] S. Eliseev, et al. Phys. Rev. Lett. 110, 082501 (2013).
[3] J. Karthein, et al. Hyperfine Interact (2019) 240: 61.

### Introduction
The following code was written to reconstruct raw PI-ICR data, to fit PI-ICR position information and calculate a frequency using the patter 1/2 scheme described in Ref. 1 and to determine a frequency ratio between a measurement ion and a reference ion. Additionally, the code allows to analyze isomeric states separated in pattern 2.

### Required software and libraries
The following code was written in Python 3.7. The required libraries are listed below with a rough description for their task in the code. It doesn't claim to be a full description of the library.

- pandas (data storage and calculation)
- numpy (calculation)
- matplotlib (plotting)
- mle (maximum likelihood estimation)
- scipy (chi square fitting)
- configparser (configuration file processing)
- jupyter (Python notebook environment)

The Python-mle package can be derived from https://github.com/ibab/python-mle, all other packages can be fetched using pip:

!pip3 install pandas numpy matplotlib scipy configparser jupyter
