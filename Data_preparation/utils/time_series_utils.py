import numpy as np

def compute_ratios(ts_dict, rec_col, treated, ts_bin):
    other = {x: [y for y in treated if y != x][0] for x in treated}
    ts_other = ts_dict[other[rec_col]]
    ts_rec = ts_dict[rec_col]
    load_rec = ts_rec[ts_bin]
    load_other = ts_other[ts_bin]
    load_tot = load_rec + load_other
    if load_tot == 0:
        return None,(None,None)
    ratio = load_rec / load_tot
    diff = np.abs(load_rec - load_other)
    winner = 'highest' if load_rec >= load_other else 'lowest'

    return ratio, (diff, winner)

class time_series_loader:
    def __init__(self, path, sep=","):
        self.ts_dict = dict()
        self.ant_dict = dict()
        with open(path) as fin:
            for row in fin:
                rows = row.strip().split(sep)
                dish_tr = rows[0]
                dish_rep = int(rows[1])
                ant_col = rows[2]
                ant_spore = rows[3]
                ant_dose = rows[4]
                ts = [float(x) for x in rows[9:]]
                ts = [0 if np.isnan(x) else x for x in ts]

                ant_tuple = (dish_tr, dish_rep, ant_col)
                if ant_spore == 'N':
                    continue
                self.ts_dict[ant_tuple] = ts
                self.ant_dict = {'spore': ant_spore,
                                 'dose': ant_dose}

    def get_ts(self, dishtr, dishrep, antcol):
        ant_tuple = (dishtr, dishrep, antcol)
        return self.ts_dict.get(ant_tuple, None)

    def get_ant_tuples(self):
        return list(self.ts_dict.keys())

    def get_ant_infto(self, dishtr, dishrep, antcol):
        ant_tuple = (dishtr, dishrep, antcol)
        return self.ant_dict.get(ant_tuple, None)

    def __getitem__(self, item):
        return self.get_ts(*item)
