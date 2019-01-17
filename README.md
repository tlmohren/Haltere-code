# Haltere project code
MATLAB code for the paper "Coriolis and centrifugal forces drivehaltere deformations and influence spike timing", by T.L. Mohren, T.L. Daniel, A.L. Eberle, P.G. Reihall, and J.L.Fox (in preparation).
This code analyzes data from FEA simulations

# Organization:
The code and data is organized through four main folders:

	figs
	functions
	scripts
	data

The scripts folder contains all the scripts that were used to create the figures and saves the figures in the figs folder.
The data is accessible through Zenodo (data found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2542944.svg)](https://doi.org/10.5281/zenodo.2542944)). If you save the data in the data folder, the scripts will be able to access it. The function addPathFolderStructure.m is called by every script to ensure access to the other folders.

# Requirements
The code was created and tested on MATLAB 2017b
