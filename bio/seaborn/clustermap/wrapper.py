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

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"]))) != "":
    makedirs(outdir)
    logging.debug(f"Directory: '{outdir}' created.")

conditions = snakemake.params["conditions"]

# Load normalized counts
data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    header=0,
    index_col=0
)

# Remove possible text annotations and validate
data = data[list(data.select_dtypes(include=[numpy.number]).columns.values)]
logging.debug("Loaded dataset:")
logging.debug(data.head())


if (nbs := len(data.columns.tolist())) != (nbc := len(list(condition.keys()))):
    message = (
        f"Expected same number of samples and conditions, got {nbs} != {nbc}"
    )
    logging.error(message)
    raise ValueError(message)

# Create custom colormap for heatmap values
cmap = seaborn.diverging_palette(
    h_neg=240,
    h_pos=10,
    as_cmap=True
)

# Create a categorical palette for samples identification
condition_lut = dict(
    zip(
        map(str, set(conditions.values())),
        seaborn.husl_palette(len(set(condition.values())), s=0.45)
    )
)


# Sorry for that part, yet I could not manage to find any
# other way to perform quicker multi-level indexing
data = data.T
data["Conditions"] = [conditions[idx] for idx in data.index]
data = (data.reset_index()
            .set_index(["Conditions", "index"])
            .T)

# Convert the palette into vectors
condition_colors = (pandas.Series(conditions.values(), index=conditions.keys())
                          .map(condition_lut))

data = data.corr()
logging.debug("Correlation table:")
logging.debug(data.head())
logging.debug("Conditional colors:")
logging.debug(condition_colors)

# Build graph
ax = seaborn.clustermap(
    data,
    cmap=cmap,
    row_colors=(
        condition_colors
        if snakemake.params.get("row_condition_color", True) is True
        else None
    ),
    col_colors=(
        condition_colors
        if snakemake.params.get("col_condition_color", True) is True
        else None
    ),
    row_cluster=(snakemake.params.get("row_cluster", True) is True),
    col_cluster=(snakemake.params.get("col_cluster", True) is True),
    linewidths=0.5,
    figsize=(
        (15 if len(data.columns.tolist()) <= 50 else 30),
        (15 if len(data.index.tolist()) < 50 else 30)
    ),
    robust=(snakemake.params.get("robust", True) is True)
)

# Rotate sample id to make them readable
matplotlib.pyplot.setp(
    ax.ax_heatmap.yaxis.get_majorticklabels(),
    rotation=snakemake.params.get("ylabel_rotation", 90)
)

matplotlib.pyplot.setp(
    ax.ax_heatmap.xaxis.get_majorticklabels(),
    rotation=snakemake.params.get("xlabel_rotation", 0)
)

# Save result
matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.info("Process over")
