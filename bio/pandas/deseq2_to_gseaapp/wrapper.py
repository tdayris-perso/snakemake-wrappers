#!/usr/bin/python3.8
# conding: utf-8

"""
Filter a DESeq2 output file
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import pandas
import numpy

def get_fc_cluster(value: numpy.float,
                   threshold: numpy.float = numpy.float(0.01)) -> str:
    """
    This functon returns a class for a given log2(Fold Change)
    """
    if abs(value) < threshold:
        return "Non-Significative"
    if value > 0:
        return "Up-regulated"
    return "Down-regulated"


def get_alpha_cluster(value: numpy.float,
                      threshold: numpy.float = numpy.float(0.05)) -> str:
    """
    This function returns a class for a given adjusted pval
    """
    if value < threshold:
        return "Differentially Expressed"
    return "Non-Significative"

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

data = pandas.read_csv(
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=0,
    dtype={
        1: str,
        2: numpy.float,
        3: numpy.float,
        4: numpy.float,
        5: numpy.float,
        6: numpy.float,
        7: numpy.float,
    }
)

padjthreshold = float(snakemake.params.get("alpha", 0.05))
fc_threshold = float(snakemake.params.get("fold_change", 0.01))
general_table_only = snakemake.params.get("general_table_only", False)

data = data[["log2FoldChange", "padj"]]

data["Cluster_FC"] = [
    get_fc_cluster(fc, fc_threshold)
    for fc in data["log2FoldChange"]
]

data["Cluster_Sig"] = [
    get_alpha_cluster(padj, padjthreshold)
    for padj in data["padj"]
]

if "fc_sig" in snakemake["output"].keys():
    logging.debug("Prining the log2(FC) / Significance table")
    tmp = data[["log2FoldChange", "Cluster_Sig"]]
    tmp.reset_index(inplace=True)
    tmp.columns = ["gene_name", "stat_change", "cluster"]
    tmp.to_csv(
        snakemake.output.fc_sig,
        sep="\t",
        index=False
    )

if "fc_fc" in snakemake["output"].keys():
    logging.debug("Prining the log2(FC) / FC cluster table")
    tmp = data[data["Cluster_Sig"] != "Non-Significative"]
    tmp = tmp[["log2FoldChange", "Cluster_FC"]]
    tmp = tmp[tmp["Cluster_FC"] != "Non-Significative"]
    tmp.reset_index(inplace=True)
    tmp.columns = ["gene_name", "stat_change", "cluster"]
    tmp.to_csv(
        snakemake.output.fc_fc,
        sep="\t",
        index=False
    )

if "padj_sig" in snakemake["output"].keys():
    logging.debug("Prining the adjusted P-Value / Significance table")
    tmp = data[["padj", "Cluster_Sig"]]
    tmp.reset_index(inplace=True)
    tmp.columns = ["gene_name", "stat_change", "cluster"]
    tmp.to_csv(
        snakemake.output.padj_sig,
        sep="\t",
        index=False
    )

if "padj_fc" in snakemake["output"].keys():
    logging.debug("Prining the adjusted P-Value / FoldChange table")
    tmp = data[data["Cluster_FC"] != "Non-Significative"]
    tmp = tmp[["padj", "Cluster_FC"]]
    tmp.reset_index(inplace=True)
    tmp.columns = ["gene_name", "stat_change", "cluster"]
    tmp.to_csv(
        snakemake.output.padj_fc,
        sep="\t",
        index=False
    )

if "complete" in snakemake["output"].keys():
    logging.debug("Prining the complete table")
    tmp = data[["log2FoldChange", "padj", "Cluster_FC", "Cluster_Sig"]]
    tmp.reset_index(inplace=True)
    tmp.columns = ["gene_name", "stat_change", "padj", "cluster", "significance"]
    tmp.to_csv(
        snakemake.output.complete,
        sep="\t",
        index=False
    )
