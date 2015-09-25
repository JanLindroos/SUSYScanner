# SUSYScanner: Python wrapper for scanning supersymmetric parameter spaces

This is a quick step quide for setting up the SUSY model and running simple 
scans over different parameter spaces

#INSTALL SUSY MODEL
In order to scan over SUSY models, the neccessary tools needs to be installed first. 
To do this go into models/SUSY and run the installation script. At present it the 
following is required before running the script:
	
	- ROOT is installed on the system with it's python module and has been set up (source thisroot.sh). 
	- C++ and fortran compilers are installed (only tested with gcc and gfortran). 
	- python is installed with modules scipy and multiprocessing
	
to run the default installation do:

	./install make

this will install softsusy, susyhit, higgsbounds (with feynhiggs), micromegas, pythia (with HEPMC), delphes and prospino
Extaernal versions of feynhiggs/HEPMC can be used for higgsbounds/pythia by specifying appropriate flags, and if some
installations fail this can be done manually by doing ./install extract [toolname] which will extract the tool in the correct 
location for manual install. Run ./install help for more details

After installation is finished everything should be set up for running a SUSY scan

#RUNNING THE CODE
All scans are run through the use of scan_cards by running run_scan.py from the main directory

	python run_scan.py [path to scan_card]
	
Choice of algorithm, number of chains, model, event generation details, etc are all set from the steering file. Some examples can be found in 
the examples folder. The scan card should be placed in the folder where you want the data from the scan to be placed. For this release only the
rand and nonrand algorithms have been tested properly, using multiprocessing paralellism (Although MPI should also work if set up with pyMPI on the cluster)

The main options to be set in the steering file are:
	
	model: This should be the path to the model file you want to run over. 4 SUSY models have been implemented (pMSSM,cMSSM, mGMSB and AMSB) with modelfiles located in models/SUSY
	alg: This is the choice of algorithm. The current options are nonrand (which takes input models from a dat file or arrays), rand (which does a random scan in the range specified in scan_range)
	mode: multiprocessing or MPI, determins wether to use shared or distributed memory pararlellization
	chains: number of parallel processes
	
In addition one can make model specific changes through model_change['model_attribute'], The most important for SUSY are:
		model_change['parts']=['softsusy','susyhit',...] specifies which tools should be run, default is all tools
		model_change['constraints']=constraints object imported from python file that contains the experimental constraints. Must be imported (see examples/pmssm_nonrand for example on how these are imported) 

	
If and when you run into bugs, please report this at https://github.com/JanLindroos/SUSYScanner or by mail to jan.lindroos@ift.uib.no.

GOOD LUCK!



