#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a clustered heatmap of multiple sample
"""

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from os.path import basename, dirname
from snakemake.utils import makedirs

# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"]))) != "":
    makedirs(outdir)

data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    index_col=0,
    header=0
)

seaborn.set(
    style="ticks",
    color_codes=True
)

g = seaborn.pairplot(
    data,
    diag_kind="kda",
    diag_kws=dict(shade=True)
)

matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
