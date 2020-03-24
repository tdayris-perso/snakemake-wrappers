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

conditions = pandas.DataFrame.from_dict(
    snakemake.params["conditions"],
    orient="index"
)
cond_id = snakemake.params.get("factor", "Condition")
conditions.columns = [cond_id]

# Load normalized counts
data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    header=0,
    index_col=0
)

# Remove possible text annotations and validate
data = data[list(data.select_dtypes(include=[numpy.number]).columns.values)]
nbs = len(data.columns.tolist())
nbc = len(conditions.index.tolist())
if nbs != nbc:
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
cond_set = set(conditions[cond_id])
colors = seaborn.husl_palette(len(cond_set), s=0.45)
cond_colors = {
    str(cond): color for cond, color in zip(list(cond_set), list(colors))
}
sample_colors = {
    sample: cond_colors[cond]
    for sample, cond in zip(conditions.index, conditions[cond_id])
}


# Sorry for that part, yet I could not manage to find any
# other way to perform quicker multi-level indexing
data = data.T
data = pandas.merge(
    data,
    conditions,
    left_index=True,
    right_index=True,
    how="left"
)
data = (data.reset_index()
            .set_index([cond_id, "index"])
            .T)

condition_colors = (pandas.Series(data.columns.get_level_values(cond_id),
                                  index=data.columns)
                          .map(cond_colors))

data = data.corr()
logging.debug("Correlation table:")
logging.debug(data.head())

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
