# collective_care

## Overview
Scripts and sample data to: 

1. Infer from data and simulate collective sanitary care of ants in the presence of a pathogen
2. Simulate a simpler exploration-exploitation model for various group sizes
 
by Casillas-Pérez B, Boďová K, Grasse AV, Tkačik G, Cremer S
Dynamic pathogen detection and social feedback shape collective hygiene in ants


## Part 1: Model inference and simulation of collective sanitary care
The subfolder ModelInference_Simulation/data contains the sample data from the collective sanitary care experiment with 6 ants. Two of the ants were treated with high fungal load (denoted by 'H' or 'F'), low fungal load ('L' or 'f'), or they were control-treated ('C' or 'Tx'), the remaining four ants were untreated. The sample data contain a single annotated replicate experiment of duration 90 minutes of each treatment group: FF, Ff, ff, FC, fC, CC.  The data contain annotated behavioral states including allogrooming and selfgrooming. The file initial_load_MM_30s.csv contains back-computed initial spore loads of all spore-treated ants using a Michaelis-Menten dynamics. The Michaelis-Menten model parameters were fitted such that the data are consistent with empirical data, as presented in the manuscript.

Note: The individual and group treatments were renamed in the process of developing the code. Some parts of the code in this repository reflect this change. 'F' for high fungal load, in the manuscript correspond to 'H', 'f' (small caps) to 'L', 'C' for control-treated to 'Tx'. The experimental group name is still constructed by combining the letters for the individual treatment (e.g.'HL' changed to 'Ff').

The subfolder ModelInference_Simulation/code contains the main script main_code.m executable in MATLAB_R2016b or later versions with a statistical toolbox. The software has been tested on OS X Yosemite, Version 10.10.5 and on Windows 10. To run the main file library minFunc_2012 and the subfolders need to be added to the Matlab path. If necessary, change the paths in the m-file load_libraries, which import these libraries. Locate the file main_code and open it in MATLAB_R2016b or newer versions. Before running the code, select the variable Fpath (line 15) as the path where the folder ModelsInference_Simulation is located. 

### System requirements
Requirement: MATLAB_R2016b with a statistical toolbox
Library minFunc_2012 for the model inference
No particular dependencies on software or operating system
The software has been tested on OS X Yosemite, Version 10.10.5 and on Windows 10
The software has been tested with MATLAB_R2017a.
No required non-standard hardware


### Installation guide
Instructions
STEP 1:
Add minFunc_2012 and the subfolders to the Matlab path. This library is used for the model inference. The full library can be found at: https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html. If necessary, change the paths in the m-file load_libraries, which import the libraries needed.

STEP 2:
Locate the file main_code and open it in MATLAB_R2016b or newer versions. Before running the code, select the variable Fpath (line 15) as the path where the folder ModelInference_Simulation is located. 

STEP 3:
Run main_code in the folder ModelInference_Simulation/code folder in MATLAB_R2016b or newer versions.

Typical install time on a "normal" desktop computer 5 min.


### Instructions to run on data

Note: In the data files, we refer to the different individual treatments by 'H' for the high spore load treatment, 'L' for the low load, and 'Tx' for the control treatment (of TritonX only). In the corresponding manuscript, we use 'F' for high fungal load (equivalent to 'H' in the data), 'f' (small caps) for low fungal load (equivalent to 'L' in the data) and 'C' for control (equivalent to 'Tx’).

Run main_code to perform three tasks:
1. Reads and statistically analyzes the exemplary input files, consisting of one experimental replicate per treatment (FF, Ff, FC, ff, fC, CC).
2. Stores the outcome of the analysis in a form of sufficient statistics in the folder output.
2. Uses the output of the statistical analysis to infer the rates of the optimal model. There are multiple classes of models that can be inferred using the code, the default is set to the RL model. The rates are visualized in a figure.
3. Uses the inferred rates from the previous step along with the initial segment of the experimental data to initialize and run a stochastic simulation of the best model (RL). The outcome of the simulation is visualized in a figure.

Parameters of the code:
line 16: run_get_stats
set to 1 to run the statistical analysis or 0 not to run it. The outcome is stored as a mat file in a folder SI_code/output/DT_statistics. The code plots the ethograms for all treatment/replicate combinations (allogrooming-green, selfgrooming-yellow). 

Optional: (lines 23, 24) the memory parameters DTA (of the P or R variable, depending on the model) and DTL (of the load variable) can be changed. Note that the statistical analysis needs to be run before the inference and the simulation with the same set of the memory parameters (which use the units 1/15 second). The memory parameters used are the ones optimal for the RL model.

line 17: run_infer_model
set to 1 to run the inference or 0 not to run it. This part of the code is quick <1s. The default model is the RL model (the optimal model in our analysis). The outcome of this part are the rates, plotted in one of the figures. The sufficient statistics (number of transitions and time spent in each bin), computed in part 1, are plotted as well, at least for the constant model and the 1D models. Note that for a better resolution one needs more replicates.

Optional: To change the model, set the variable “variant” (line 167) to a value 1-9:
1=constant model
2=P model
3=R model
4=La model
5=Las model
6=PLa model
7=PLas model
8=RLa model
9=RLas model

line 18: run_simulate
set to 1 to run the simulation corresponding to a single replicate from the provided empirical data. The simulation runs model RL (the best model). The outcome of this part is composed of two ethograms: one experimental, one simulated. The simulated one will have only 30 min worth of data unless changed in the code (the simulation slows down when larger times are used).

Optional: The simulation uses the first treatment and the first replicate. To change this modify idtr in line 318 (to change the treatment, allowed values 1-6, 1=FF, 2=Ff, 3=FC, 4=ff, 5=fC, 6=CC). The variable idrep (to change the replicate) in line 319 can be changed only when more experimental replicates are loaded.

Optional: The simulation generates sanitary behavior corresponding to the first 30 min of the experiment. This interval can be changed in line 288 by setting a different value of T (T=15*60*90 corresponds to 90 min). Note that the longer time results in a longer computation time.
 
Expected runtime:
- Obtaining sufficient statistics from the data: <3 min when all 6 supplied input files are analyzed.
- Inferring transition rates: <1s
- Stochastic simulation of one replicate for 30 min: <2 min
- Stochastic simulation of one replicate for 90 min: <5 min
 
Expected output:
- sufficient statistics stored in ModelInference_Simulation/output/DT_statistics as a mat file
- inferred transition rates, visualized in a figure
- output of the stochastic simulation plotted as a raster plot along with a similar graph for the corresponding experimental data


### Additional instructions for use
The code is capable of running simple examples as well as the complex tasks in the manuscript. To obtain the key results in our manuscript one needs to supply all experimental replicates and set the parameters appropriately.

### File description for the code items


main_code.m:
Runs all three major parts of the code:

1. Loads the experimental data and collects sufficient statistics

2. Based on the statistics computes the transition rates using log-likelihood minimization

3. Uses the inferred rates to generate a stochastic simulation of the ant dish and compares it with the experimental data

### Code labels for part 1
importfile.m:
Imports initial loads of all considered ants.

readdish.m:
Reads the experimental data and stores it into convenient data formats: list of grooming events, their durations, loads, table of all grooming states with times organized by event type (allo/self) and ant.

getstate.m:
Classifies the experimentally recorded states into allo/self grooming.

detailed_states_trace.m:
Computes the necessary statistics in terms of time, state, and load.

get_dishstats.m:
Computes the detailed frame-by-frame vector of performed/received activity, own load, load seen during grooming from the list of grooming events with times.

getstats_0D.m:
Computes statistics (number of transitions between states, time spent in each state) independent of other variables.

getstats_2D.m:
Computes statistics (number of transitions between states, time spent in each state) with 2D dependence, where the first dependent variable is assumed to be an activity variable while the second one can be a load variable.

### Code labels for part 2
plotrates.m:
Plots the inferred rates, computed in part 2 of the main code.

shapeMtoV.m,
shapeVtoM.m:
Reshapes the transition rate for the purpose of the inference. One file reshapes them into the input format of the minimization routine, the other reshapes the outcome back to the original dimensions.

Inference_rho.m:
Infers all model parameters, with the transition rates and the parameter rho.

lik0D.m:
Definition of likelihood for the inference independent of other factors.

lik1D.m:
Definition of likelihood for the inference dependent on one factor.

lik2D.m:
Definition of likelihood for the inference dependent on two factors.

lik2DR.m:
Definition of likelihood for the inference dependent on two factors and parameter rho.

### Code labels for part 3
- all within the file main_code.m

### Auxiliary functions
load_libraries.m:
Loads the minFunc library for multidimensional optimization. The paths have to be set appropriately.

Col.m:
Defines the coloring schemes used in the plots


## Part 2: Exploration-exploitation model. 
The folder “EE_forward_simulation” contains the script and the separate description file explaining the use, the parameters, and the output of the script. 


### System requirements
Requirement: MATLAB_R2016b with a statistical toolbox
No particular dependencies on software or operating system
The software has been tested on OS X Yosemite, Version 10.10.5 and on Windows 10
The software has been tested with MATLAB_R2017a.
No required non-standard hardware

### Installation guide
Instructions
STEP 1:
Open the file exploration_exploitation in MATLAB_R2016b or newer versions. 

### Demo
The code runs the exploration-exploitation model with the main parameters:
N - number of ants
NL - number of exposed ants 
L0, ep - sequential rule parameters
K - max rule parameter - number of checked ants
tE, tG - time of exploration/grooming
ML - initial load
rV, rK - MM dynamics, 

The preset parameter values:
N = 30, NL = 10; L0 = 10, ep = 0.01, K = 30, tE = 10, tG = 100, rV = 1, rK = 100, ML = 100


Run main_code to perform the following tasks:
1. Run the exploration-exploitation model with the MAX rule with the main parameter K<=N (the number of ants explored before choosing the one with the highest load for grooming)
2. Run the exploration-exploitation model with the SEQ rule with the main parameters L0, ep
 
Expected runtime:
Depends on the parameters of the model (tE, tG, and mainly the terminal time T, which can be set separately for each model version (lines 43 and 54) but for the preset choice of parameters the run time is ~ 5s.
 
Expected output:
Figures of the results of the model corresponding to the two rules into figure 1 (MAX) and figure 2 (SEQ) in terms of the raster plots of performed and received grooming and the load intensity. The ants in the rasters are ordered by the value of the initial spore load from the highest to the lowest.
