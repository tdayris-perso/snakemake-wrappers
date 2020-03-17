"""Snakemake wrapper for SnpSift gwasCat"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import dirname
from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")
makedirs(dirname(snakemake.output.vcf))

shell(
    "SnpSift gwasCat"  # Tool and its subcommand
    " {extra}"  # Extra parameters
    " -db {snakemake.input.gwascat}"  # Path to annotation vcf file
    " {snakemake.input.vcf}"  # Path to input vcf file
    " > {snakemake.output.vcf}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)