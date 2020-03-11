#!/usr/bin/python3.8
# conding: utf-8

"""
Join multiple Salmon files, rename columns and
add genes name if annotation is provided
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import logging
import pandas
import numpy

from os.path import dirname, basename
from snakemake.utils import makedirs


def read_tx2gene(path: str,
                 genes: bool = False,
                 header: bool = False) -> pandas.DataFrame:
    """
    This function reads a TSV file containing the following columns:

    1: ens_gene
    2: ens_transcript
    3: ext_gene

    And returns it as a DataFrame
    """
    t2g = pandas.read_csv(
        path,
        sep="\t",
        index_col=1,
        header=(0 if header is True else None),
        dtype=str
    )
    t2g.columns = ["Ensembl_ID", "Hugo_ID"]

    if genes is True:
        t2g = (t2g.reset_index()[["Ensembl_ID", "Hugo_ID"]]
                  .drop_duplicates()
                  .set_index("Ensembl_ID"))

    return t2g


def read_salmon(path: str) -> pandas.DataFrame:
    """
    This function reads a single salmon quant.sf or quant.genes.sf
    and returns it as a pandas DataFrame
    """
    return pandas.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str,
            1: numpy.float,
            2: numpy.float,
            3: numpy.float,
            4: numpy.float
        },
        na_values=""
    )


if (outdir := dirname(snakemake.output["tsv"])) != "":
    makedirs(outdir)

merged_frame = None

for quant in snakemake.input["quant"]:
    data = read_salmon(quant)

    sample_id = basename(dirname(quant))
    if len(suffix := snakemake.params.get("suffix", "")) > 0:
        sample_id = sample_id[:-len(suffix)]
    if len(prefix := snakemake.params.get("prefix", "")) > 0:
        sample_id = sample_id[len(prefix):]

    data = data[[snakemake.params.get("column", "TPM")]]
    data.columns = [sample_id]

    try:
        merged_frame = pandas.merge(
            merged_frame,
            data,
            left_index=True,
            right_index=True
        )
    except TypeError:
        merged_frame = data

merged_frame.fillna(0)

if snakemake.params.get("gencode", False) is True:
    merged_frame = merged_frame.set_index(
        pandas.DataFrame(merged_frame.index.str.split(".", 1).tolist())[0]
    )

if snakemake.params.get("drop_null", False) is True:
    merged_frame = merged_frame.loc[~(merged_frame == 0).all(axis=1)]

if snakemake.params.get("drop_na", False) is True:
    merged_frame.dropna(axis=0, how="all", inplace=True)

if (tr2gene_path := snakemake.input.get("tx2gene", None)) is not None:
    merged_frame = pandas.merge(
        merged_frame,
        read_tx2gene(
            tr2gene_path,
            snakemake.params.get("genes", False),
            snakemake.params.get("header", False)
        ),
        left_index=True,
        right_index=True,
        how="left"
    )

merged_frame.to_csv(
    snakemake.output["tsv"],
    sep="\t",
    index=True,
    header=True,
    index_label=(
        "target_id"
        if snakemake.params.get("index_label", False) is True
        else False
    )
)