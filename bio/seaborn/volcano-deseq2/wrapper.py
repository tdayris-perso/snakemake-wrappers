#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a Volcano plot
"""

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from os.path import basename, dirname
from snakemake.utils import makedirs
from typing import Generator, List

# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"][0]))) != "":
    makedirs(outdir)


def classify(alphas: List[numpy.float],
             foldchanges: List[numpy.float],
             alpha_threshold: numpy.float,
             fc_threshold: numpy.float) -> Generator[str, None, None]:
    """
    This function classifies the targets among three conditions:
    1: Differentially Expressed
    2: Below alpha threshold
    3: Below fc threshold
    """
    for alpha, foldchange in zip(alphas, foldchanges):
        if alpha <= alpha_threshold:
            if abs(foldchange) <= fc_threshold:
                yield "Differentially Expressed"
            else:
                yield "Below fold change threshold"
        else:
            yield "Below significance threshold"


# Load dataset
data = pandas.read_csv(
    snakemake.input["deseq2_tsv"],
    sep="\t",
    index_col=0,
    header=0,
    dtype=numpy.float
)

# Filter content
data = data[["log2FoldChange", "padj"]]
data.dropna(axis=0, how="all", inplace=True)
data = data.loc[~(data == 0).all(axis=1)]

# Define filters' thresholds
alpha_threshold = numpy.log(snakemake.params.get("alpha_threshold", 0.05))
fc_threshold = -numpy.log10(snakemake.params.get("fc_threshold", 0.001))
data["Significance"] = list(
    classify(
        alphas=data.padj,
        foldchanges=data.log2FoldChange,
        alpha_threshold=alpha_threshold,
        fc_threshold=fc_threshold
    )
)

# Transform padj in order to plot it
data["-log10(padj)"] = [
    0 if padj == 0 else -numpy.log10(padj) for padj in data.padj
]

# Build plot
seaborn.set(
    style="ticks",
    palette="pastel",
    color_code=True
)

g = seaborn.lmplot(
    x="log2FoldChange",
    y="-log10(padj)",
    hue="Significance",
    hue_order=[
        "Below significance threshold",
        "Below beta threshold",
        "Differentially Expressed"
    ],
    data=data,
    palette="Blues",
    fit_reg=False
)

# Add filter lines
g.map(
    matplotlib.pyplot.axhline,
    y=alpha_threshold,
    ls=":",
    c=".5",
    markersize=0.1
)

g.map(
    matplotlib.pyplot.axvline,
    x=fc_threshold,
    ls=":",
    c=".5",
    markersize=0.1
)

g.map(
    matplotlib.pyplot.axvline,
    x=-fc_threshold,
    ls=":",
    c=".5",
    markersize=0.1
)

# Save final plot
matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
