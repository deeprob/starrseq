import os
import sys
import numpy as np
import pandas as pd
import subprocess
import itertools
from itertools import starmap
from argparse import Namespace

from bokeh.layouts import row,column
from bokeh.plotting import figure, show
from bokeh.palettes import Blues9, Spectral10
from bokeh.models import ColorBar, LinearColorMapper, ColumnDataSource, FactorRange, FixedTicker, NumeralTickFormatter

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

def get_tsv(filename):
    return pd.read_csv(filename, header=None, sep="\t")

def get_rep_cov_names(lib_pre, lib_rep, lib_suff):
    """Helper function to get the replicate coverage and depth files names """
    lib_reps = lib_rep.split()
    fnames = []
    for rep in lib_reps:
        fname = "_".join([lib_pre, rep, lib_suff])
        fnames.append(fname)
    return fnames

def get_corr_matrix(dfs, col_idx=3):
    """Correaltion matrix between replicates"""
    counts_col = col_idx
    mat = np.ones((len(dfs), len(dfs)))
    mat[0,1] = mat[1,0] = dfs[0][counts_col].corr(dfs[1][counts_col])
    mat[0,2] = mat[2,0] = dfs[0][counts_col].corr(dfs[2][counts_col])
    mat[1,2] = mat[2,1] = dfs[1][counts_col].corr(dfs[2][counts_col])
    return mat

def rep_bed_to_df(lib_aligned_prefix, lib_reps, lib_suffix):
    """Replicate coverage bed file to dataframe objects"""
    rep_cov_beds = get_rep_cov_names(lib_aligned_prefix, lib_reps, lib_suffix)
    rep_cov_dfs = list(map(get_tsv, rep_cov_beds))
    return rep_cov_dfs

def df_to_corr_mat(dfs, col_idx):
    """Dataframe objects to correlation matrix"""
    return get_corr_matrix(dfs, col_idx)

def rep_bed_to_corr_mat(lib_aligned_prefix, lib_reps, lib_suffix, col_idx):
    rep_cov_dfs = rep_bed_to_df(lib_aligned_prefix, lib_reps, lib_suffix)
    return df_to_corr_mat(rep_cov_dfs, col_idx)

def round_n(val):
       return round(val, 3)

def get_intra_lib_rep_corr_plot(in_corrmat, cc_corrmat, ko_corrmat, fig_title="Intra library replicate correlation of read depths"):
    arr_in = in_corrmat.flatten(order="C")
    arr_cc = cc_corrmat.flatten(order="C")
    arr_ko = ko_corrmat.flatten(order="C")
    mapper = LinearColorMapper(palette=Blues9, low=min(arr_in.min(), arr_cc.min(), arr_ko.min()), high=1)

    TOOLS = "save,pan,box_zoom,reset,wheel_zoom"

    exp = ["Input", "Control", "KO"]
    reps = ["Rep_1", "Rep_2", "Rep_3"]
    exp_reps = list(itertools.product(exp, reps))

    mydict = {"a":sum([list(itertools.repeat(a, 3)) for a  in exp_reps], []),
              "b":sum(itertools.repeat(reps, 9),[]),
              "values":list(map(round_n, arr_in)) + list(map(round_n, arr_cc)) + list(map(round_n, arr_ko))}

    source = ColumnDataSource(mydict)


    p = figure(title=fig_title,
               x_range=FactorRange(*exp_reps), y_range=FactorRange(*reps),  
               x_axis_location="below", width=900, height=300,
               tools=TOOLS, toolbar_location='left')

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "7px"
    p.axis.major_label_standoff = 0

    p.rect(x="a", y="b", width=0.95, height=0.95, source=source,
           fill_color={"field":"values", "transform": mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
                         label_standoff=5, border_line_color=None)
    p.add_layout(color_bar, 'right')

    p.text(mydict["a"], mydict["b"], text=mydict["values"],
           text_baseline="middle", text_align="center")

    show(p)
    return



def get_fractional_base_coverage_plot(cov_df, col_idx=6, lib_name="input"):
    arr = np.histogram(cov_df[col_idx])
    TOOLS = "save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title=f"Fraction of the genomic regions of interest covered \nby the {lib_name} library", 
               x_axis_location="below", width=400, height=400,
               tools=TOOLS, toolbar_location='left')

    mydict = dict(fraction=arr[1][:-1] + 0.05, counts=arr[0], color=Spectral10)
    source = ColumnDataSource(data=mydict)

    p.grid.grid_line_color = None
    p.axis.major_label_text_font_size = "6px"
    p.axis.major_label_standoff = 0
    p.xaxis.ticker = FixedTicker(ticks=arr[1])

    p.vbar(x="fraction", top="counts", width=0.095, source=source, color="color")

    p.text(mydict["fraction"], mydict["counts"] + 1000, text=mydict["counts"],
           text_baseline="middle", text_align="center")

    return p

def get_inter_lib_scatter_plot(in_df, cc_df, ko_df, col_idx=3):
    s1 = figure(title=f"Control vs Input library read counts across regions of interest \nPearson Correlation: {round_n(in_df[3].corr(cc_df[3]))}",
                width=450, height=350, background_fill_color="#fafafa")
    s1.circle(in_df[col_idx], y=cc_df[col_idx], size=12, color="firebrick", alpha=0.8)
    s1.xaxis.axis_label = 'Input library read counts'
    s1.yaxis.axis_label = 'Control library read counts'

    s2 = figure(title=f"KO vs Control library read counts across regions of interest \nPearson Correlation: {round_n(cc_df[3].corr(ko_df[3]))}",
                width=450, height=350, background_fill_color="#fafafa")
    s2.circle(cc_df[col_idx], y=ko_df[col_idx], size=12, color="darkgreen", alpha=0.8)
    s2.xaxis.axis_label = 'Control library read counts'
    s2.yaxis.axis_label = 'KO library read counts'

    s3 = figure(title=f"KO vs Input library read counts across regions of interest \nPearson Correlation: {round_n(in_df[3].corr(ko_df[3]))}",
                width=450, height=350, background_fill_color="#fafafa")
    s3.circle(in_df[col_idx], y=ko_df[col_idx], size=12, color="navy", alpha=0.8)
    s3.xaxis.axis_label = 'Input library read counts'
    s3.yaxis.axis_label = 'KO library read counts'

    show(row(s1, s2, s3))   
    return

def get_hist_info(df):
    return np.histogram(df[4], bins=[0, 1, 10, 20, 40, 50, 100, 1000, 10000, 100000, 1000000])

def get_bps_covered(hist_info, readcount_val):
    return sum(hist_info[0][np.argwhere(hist_info[1]>=readcount_val)[0][0]:])

def get_lib_percent_covered(df, readcount_val=50):
    hist_info = get_hist_info(df)
    bps_covered = get_bps_covered(hist_info, readcount_val)
    return bps_covered/len(df), bps_covered, len(df)

def print_library_bp_coverage_stats(df, lib_name, readcount_val=50):
    lib_stats = get_lib_percent_covered(df, readcount_val)
    print(f"----{lib_name} library stats----\n")
    print(f"Number of Bases which have at least {readcount_val} reads assigned to it: {lib_stats[1]}")
    print(f"Total number of Bases in the library: {lib_stats[2]}")
    print(f"Percentage of Bases which have greater than {readcount_val} reads assigned to it: {lib_stats[0]}\n")
    return

def get_cradle_activated(cradle_peak_file, cradle_activated_file, effect_size=0):
    cradle_df = pd.read_csv(cradle_peak_file, sep="\t")
    cradle_df.loc[cradle_df["effectSize"]>effect_size].iloc[:, [0, 1, 2, 3, 4, 5, 6, 7]].to_csv(cradle_activated_file, header=None, index=None, sep="\t")
    return 

def get_cradle_repressed(cradle_peak_file, cradle_repressed_file, effect_size=0):
    cradle_df = pd.read_csv(cradle_peak_file, sep="\t")
    cradle_df.loc[cradle_df["effectSize"]<effect_size].iloc[:, [0, 1, 2, 3, 4, 5, 6, 7]].to_csv(cradle_repressed_file, header=None, index=None, sep="\t")
    return 

def get_mea(activity_file, reference_genome, background_region, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    logfile = os.path.join(output_dir, "homer.log")

    lf = open(logfile, "w")
    subprocess.run(["bash", f"{CURRENT_DIR_PATH}/run_mea.sh", activity_file, 
                    reference_genome, background_region, output_dir], stdout=lf, stderr=lf, check=True)
    lf.close()
    return logfile
