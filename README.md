# CMIP6 Arctic Terns (replace with actual title when we choose it...)
This is the github repository for scripts required to reproduce figures and analyses in Morten et al. (202X) [Title].

Repository structure:

cmip6_arctic_terns/
├─ ternsim/                      
│  ├─ PREPROCESSING.sh              | Shell script used to preprocess raw CMIP6 output for particle tracking
│  ├─ TernSim.py                    | Simulates vTern trajectories from preprocessed CMIP6 surface velocity data
│  ├─ flight_path.py 	              | Generates Figure 4
│  ├─ flight_time.py 	              | Generates Supplementary Figure 3
│  ├─ path_comp.py	                | Generates Supplementary Figure 2
│  ├─ ternmethods.py                | Methods for TernSim.py   
│  ├─ FIGURES/
│  │  ├─ vTern_path_future.pdf      | Figure 4
│  │  ├─ vTern_path_validation.pdf  | Supp. Figure 2
│  │  ├─ vTern_time_future.pdf      | Supp. Figure 3
│  ├─ resources/
│  │  ├─ gridsource.txt             | Regridding file used by PREPROCESSING.sh



