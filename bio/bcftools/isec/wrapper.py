"""Snakemake wrapper for BCFTools isec"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from os.path import dirname
from snakemake.utils import makedirs
from snakemake.shell import shell

makedirs(dirname(snakemake.output[0]))

shell(
    "bcftools isec "  # BCFTools and its subcommand
    "{snakemake.params} "  # Extra parameters
    "-o {snakemake.output[0]} "  # Path to output file
    "{snakemake.input.calls}"  # Path to input VCF/BCF files
)
