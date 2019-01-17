# Haltere project code
	Code for the paper "Coriolis and centrifugal forces drive haltere deformations and influence spike timing", by
	T.L. Mohren, T.L. Daniel, A.L. Eberle, P.G. Reihall, and J.L. Fox (in preparation).
	This code analyzes data from FEA simulations (see xxxx).

# Organization:
The code and data is organized through four main folders:

	figs
	functions
	scripts
	data

The scripts folder contains all the scripts that were used to create the figures and saves the figures in the figs folder.
When cloning this repository, you will not have the data folder.
The data is accessible through Zenodo (xxxx). If you create a data folder alongside the other folders and save the data in there, the scripts will be able to access it. The function addPathFolderStructure is called by every script to ensure access to the other folders. 
