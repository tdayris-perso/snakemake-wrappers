__author__ = "Jan Forster"
__copyright__ = "Copyright 2019, Jan Forster"
__email__ = "jan.forster@uk-essen.de"
__license__ = "MIT"

import os

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell("gatk --java-options '{java_opts}' SplitNCigarReads {extra} "
    " -R {snakemake.input.ref} -I {snakemake.input.bam} "
    "-O {snakemake.output} {log}")
