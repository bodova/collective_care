
# Create analysis ready tables from solomon and ddpcr data
# Usage: execfile('sampleAnalysis.py') or python scripts/sampleAnalysis.py

import sys
import os
sys.path.append(os.getcwd())

import matplotlib
from matplotlib import pylab as P

import utils.tableBuilder as tB
from utils.analysis_functions import plotEventDurations
from utils.choice_distribution import choice_distribution_plots
import configs.config_tests as conf

matplotlib.rcParams.update({'font.size': 22})

# Initialize object
if any(["table" in x for x in conf.steps_to_run]):
    print("Creating table reader object.")
    mybuilder = tB.TableBuilder(solomonDir=conf.solomonDir,
                                outFolder=conf.tableOutput,
                                ddpcrFilePath=conf.ddpcrFilePath,
                                metadataFilePath=conf.metadataFilePath)

if "windowed_tables" in conf.steps_to_run:
    print("Building windowed tables")
    for beh in tB.behaviourTypes:
        mybuilder.buildWindowedTables(treatmentReplicatePairs=conf.treatmentReplicatePairs,
                                  windowSizeInFrames=conf.windowSizeInFrames,
                                  behaviourType=beh,
                                  eventStartCount=conf.eventStartCount)

if "windowed_tables_pairwise" in conf.steps_to_run:
    print("Building pairwise windowed tables")
    mybuilder.buildWindowedTablesPairwise(conf.treatmentReplicatePairs,
                                windowSizeInFrames=conf.windowSizeInFrames,
                                eventStartCount=conf.eventStartCount)

if "build_all_tables" in conf.steps_to_run:
    print("Building all tables")
    mybuilder.buildAllTables(conf.treatmentReplicatePairs)


if "plot_choice_distribution" in conf.steps_to_run:
    print("Plotting histograms of choice")
    choice_distribution_plots(windowSizeInFrames=conf.windowSizeInFrames,
                              MM_files_location=conf.MM_data_path,
                              treatments=conf.treatments,
                              data_files_location=conf.tableOutput,
                              output_files_location=conf.tableOutput
                              )




