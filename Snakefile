import pandas as pd

configfile: "config.yaml"
working_dir: config["outdir"]
samples_df = pd.read_csv(config["samples"], sep="\t", comment="#")

include: "rules/download.smk"
include: "rules/qc.smk"


rule all:
    input:
        rules.multiqc.output