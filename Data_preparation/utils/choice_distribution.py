import os
import random as ran
import glob

import numpy as np
from matplotlib import pylab as P

from utils.time_series_utils import time_series_loader
from utils.time_series_utils import compute_ratios


def choice_distribution_plots(windowSizeInFrames,
                              treatments,
                              data_files_location,
                              MM_files_location,
                              output_files_location,
                              max_frame=85 * 60 * 15,
                              num_nulls=100,
                              numbinsx=20,
                              numbinsy=20,
                              uselog_in2d=False,
                              normalize_loads=True,
                              maxeventduration=np.inf,
                              show_low_choices=True,
                              plot_vs_null=False,
                              plot_diff_vs_duration=False
                              ):
    use_ratios_for2d = False
    use_durations_for_bars = False
    use_log_for_bars = True

    figs = []

    MM_loads_file = os.path.join(MM_files_location, "POST_windowed_" + str(
        windowSizeInFrames) + "_groomTot_spore_MM_paper.csv")
    eventfiles = glob.glob(data_files_location + "EVENTS_POST_*")
    eventfiles.sort()

    all_ts = time_series_loader(MM_loads_file)
    num_windows = len(all_ts[all_ts.get_ant_tuples()[0]])

    window_edges = [i * windowSizeInFrames for i in range(num_windows + 1)]
    ant_colors = set([x[2] for x in all_ts.get_ant_tuples()])

    event_load_ratios = []
    event_loads = []
    event_durations = []
    null_load_ratios = []
    null_loads = []
    events_with_differences = []

    #  0  1      2     3      4     5  6  7     8    9  10
    # TxTx,  1, 21978, 47,   yellow, T, X, S,  yellow, T, X
    # dtr, rep   t    dur    acol   ,  , typ,  rcol , ,
    fout_ratios_real = open(os.path.join(output_files_location, "ratios_real.csv"), "w")
    fout_ratios_null = open(os.path.join(output_files_location, "ratios_null.csv"), "w")
    real_headers = ["treatment", "replicate","window", "load_ratio","duration","load_diff", "target"]
    rand_headers = ["treatment", "replicate", "window", "load_ratio"]
    fout_ratios_real.write("#"+",".join(real_headers)+"\n")
    fout_ratios_null.write("#" + ",".join(rand_headers) + "\n")
    for fn in eventfiles:
        if all([prefix not in fn for prefix in treatments]):
            continue
        num_events = 0
        print(fn)
        with open(fn) as fin:
            dishtr = fn.split("/")[-1].split('.')[0].split("_")[2]
            dishrep = int(fn.split("/")[-1].split('.')[0].split("_")[3])
            these_ts = {col: all_ts[(dishtr, dishrep, col)] for col in ant_colors}
            treated = [ac for ac in ant_colors if these_ts[ac] is not None]

            for event_row in fin:
                rows = event_row.strip().split(",")
                type_event = rows[7]
                treat_r = rows[9]
                if type_event != "G" or treat_r == 'N':
                    continue

                time = int(rows[2])
                duration = int(rows[3])
                if time > max_frame:
                    continue

                # In which window is this event
                ts_bin = np.argmin(np.array([np.abs(time - bedge)
                                             for bedge in window_edges
                                             if time - bedge >= 0]))
                rec_col = rows[8]
                # these_ts: time series for this dish
                # rec_col:  who is recieving
                # treated:  all who are treated
                # ts_bin:   in which window of the timeseries is this happening
                # ratio:
                # diff_winner is a tuple (diff, winner) where diff is the abs of the difference between the two
                #                                         and winner is either 'highest' or 'lowest'
                ratio, diff_winner = compute_ratios(these_ts, rec_col, treated, ts_bin)
                if ratio is None or np.isnan((ratio)):
                    continue
                num_events += 1
                event_load_ratios.append(ratio)
                event_loads.append(these_ts[rec_col][ts_bin])
                event_durations.append(duration)

                # here we generate "random events"
                for i in range(num_nulls):
                    # for null_col in treated:
                    null_col = ran.sample(treated, 1)[0]
                    nullrat, _ = compute_ratios(these_ts, null_col, treated, ts_bin)
                    if nullrat is None or np.isnan((nullrat)):
                        i += 1
                        continue
                    null_load_ratios.append(nullrat)
                    null_loads.append(these_ts[null_col][ts_bin])
                    fout_ratios_null.write(",".join([str(x) for x in [dishtr, dishrep, ts_bin, nullrat]]) + "\n")

                diff = diff_winner[0]
                winner = diff_winner[1]
                fout_ratios_real.write(",".join([str(x) for x in [dishtr,
                                                                  dishrep,
                                                                  ts_bin,
                                                                  ratio,
                                                                  duration,
                                                                  diff,
                                                                  winner]]) + "\n")
                events_with_differences.append((diff, winner, duration, dishtr, ratio))
                # The 2D Plot is the histogram of pairs of the form (duration, load)
                #    where  load = ratio  if  normalize_odds
                #       or  load = diff if    not normalize odds
                # duration is the 5th column of fout_ratios_real file
                # ratio    is the 4th column of fout_ratios_real file
                # diff     is the 6th column of fout_ratios_real file
            print(dishtr, dishrep, num_events, treated, num_events)

    print("found", len(event_load_ratios), "grooming to treated events")
    chose_highest = len([x for x in events_with_differences if x[1] == "highest"]) / len(event_load_ratios)
    print("of these,", chose_highest * 100, "%  were grooming the ant with a higher load")

    fout_ratios_real.flush()
    fout_ratios_real.close()
    fout_ratios_null.flush()
    fout_ratios_null.close()

    if plot_vs_null:
        # -----------------------------------------------------------------------
        fi = P.figure()
        P.hist((event_loads), bins=20, density=True)
        P.hist((null_loads), bins=20, density=True, alpha=0.5)
        fi.legend(["observed", "expected"])
        P.axvline(x=0.5, color='red')
        P.title("Distribution target load at start of grooming event\nfor: " +
                " ".join([x for x in treatments]))
        P.xlabel("L_target ")

        # -----------------------------------------------------------------------
        fi2 = P.figure()
        P.hist(event_load_ratios, bins=20, density=True)
        P.hist(null_load_ratios, bins=20, density=True, alpha=0.5)
        P.legend(["observed", "expected"])
        P.axvline(x=0.5, color='red')
        P.title("Distribution target load ratio at start of grooming event\nfor: " +
                " ".join([x for x in treatments]))
        P.xlabel("L_target / (L_target + L_other)")

        figs += [fi, fi2]
    # -----------------------------------------------------------------------

    # #   2D histogram of load difference vs event duration
    if plot_diff_vs_duration:
        fi3 = P.figure()
        selectors = ["highest", "lowest"] if show_low_choices else [None]

        for nse, selector in enumerate(selectors):
            if len(selectors) > 1:
                P.subplot(1, 2, nse + 1)
            durations = [x[2] if x[0] > 0 else 0 for x in events_with_differences
                        if selector is None or x[1] == selector]
            if normalize_loads:
                loadcoord = 4
                ylab = "L_taget / (L_taget + L_other)"
            else:
                loadcoord = 0
                ylab = "|L_taget - L_other|"

            if uselog_in2d:
                loads = [np.log10(x[loadcoord]) if x[loadcoord] > 0 else 0 for x in events_with_differences
                        if selector is None or x[1] == selector]
                ylab = "log(" + ylab + ")"
            else:
                loads = [x[loadcoord] for x in events_with_differences
                        if selector is None or x[1] == selector]

            loads = [loads[i] for i, d in enumerate(durations) if d < maxeventduration]
            durations = [d for i, d in enumerate(durations) if d < maxeventduration]

            ranges = np.zeros((2, 2))
            ranges[0, 0] = min(durations)
            ranges[0, 1] = max(durations)
            ranges[1, 0] = 0 if normalize_loads else min(loads)
            ranges[1, 1] = 1 if normalize_loads else max(loads)

            P.hist2d(x=np.array(durations), y=loads,
                    bins=(numbinsx, numbinsy),
                    range=ranges)
            P.xlabel("Duration of event ")

            if len(selectors) > 1:
                P.title(selector)
            else:
                P.title("Event durations Vs load ratios at start of grooming event\nfor " + " ".join(treatments))

            if nse == 0:
                P.ylabel(ylab)
            else:
                P.yticks([])

        figs.append(fi3)


    return figs
