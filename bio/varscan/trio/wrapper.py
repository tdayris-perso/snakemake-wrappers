"""Snakemake wrapper for varscan trio"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from os.path import dirname, splitext
from snakemake.shell import shell
from snakemake.utils import makedirs

# Defining logging and gathering extra parameters
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

# Building output dirs
makedirs(dirname(snakemake.output.snp))
makedirs(dirname(snakemake.output.indel))

# Output prefix
prefix = snakemake.params.get("prefix", splitext(snakemake.output.snp)[0])

# Searching for input files
shell(
    "varscan trio"  # Tool and its subcommand
    " {snakemake.input.mpileup}"  # Path to input file(s)
    " {prefix}"  # Path to output
    " {extra}"  # Extra parameters
    " {log}"  # Logging behaviour
)
