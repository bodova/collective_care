basedir = "./data"

# Input (read-only) directories-----------------------------------

## Spore measurement dataset
ddpcrDir = basedir+"/raw/ddPCR/"
ddpcrFilePath = ddpcrDir+'/IDE_ddPCR_results_final.csv'
## Behavioral annotations dataset output by Solomon coder
solomonDir = basedir+"/raw/BehaviorSolomonData/"
## Metadata 
metadataFilePath = basedir+"/raw/BehaviorSolomonData/behavioural_coding_metadata.csv"
MM_data_path = basedir+"/intermediate/" 
# Output directories-----------------------------------------------
tableOutput = basedir+"/output/"  

# Pipeline-----------------------------------------------------
steps_to_run = [
    "windowed_tables", 
    "windowed_tables_pairwise",
    "build_all_tables",
    "plot_choice_distribution" # see NOTE 
    ] 

# NOTE: The final step in the pipeline is to run 'plot_choice_distribution',
#  which requires output generated in previous steps. Specifically, it uses files 
# produced by 'build_all_tables' (EVENTS*.csv) and a table
#  (./data/intermediate/POST_windowed_450_groomTot_spore_MM_paper.csv') containing a
#  time series of estimated current spore loads (as described in
#  Modeling_and_simulations/MM_backcomputation). The creation of this table requires
#  files generated by the 'windowed_tables' step (POST_windowed_450_groomTot.csv).



# Parameters --------------------------------------------------------
## List of replicates to be processed
treatmentReplicatePairs = [("HL", 16),
                           ("HH", 1),
                           ("TxTx", 18)]

## For windowed functions
windowSizeInFrames = 30 * 15 # sec * fps = frames

# For pair-wise timeseries
eventStartCount=False   # True :place event starts in a window, instead of the number of frames 

# For choice ratios
max_frame = 85*60*15  # Events beyond this frame are discarded
maxeventduration = 600  # Events longer than this are discarded
num_random_for_null = 1000
treatments = ["HH", "HL", "LL"]









