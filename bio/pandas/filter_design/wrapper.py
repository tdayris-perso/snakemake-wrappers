#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a design file
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import pandas
import numpy

data = (pandas.read_csv(snakemake.input["design"], sep="\t", header=0)
              .set_index("Sample_id"))

if (samples := snakemake.params.get("keep_sample", None)) is not None:
    data = data[data.index.isin(samples)]

if (not_samples := snakemake.params.get("filter_sample", None)) is not None:
    data = data[not data.index.isin(not_samples)]

if (cols := snakemake.params.get("keep_column", None)) is not None:
    data = data[cols]

if (not_cols := snakemake.params.get("filter_column")) is not None:
    data = data[list(set(data.column.tolist()) - set(not_cols))]

data.to_csv(
    snakemake.output["tsv"],
    sep="\t"
)
