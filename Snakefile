import pandas as pd
from pathlib import Path

configfile: "config.yaml"
working_dir: config["outdir"]
samples_df = pd.read_csv(config["samples"], sep="\t", comment="#")

include: "rules/download.smk"
include: "rules/qc.smk"
include: "rules/align.smk"

# Set brbseq samples
for d in config["cellbarcodes"]:
    sample = d["sample"]
    cellbarcode_file = Path(d["file"])

    assert cellbarcode_file.exists(), f"Cell barcode file {cellbarcode_file} does not exist"
    samples_df.loc[samples_df["sample"] == sample, "cellbarcode_file"] = cellbarcode_file


rule all:
    input:
        rules.multiqc.output,
        expand(
            rules.star_solo.output,
            sample=samples_df.loc[samples_df["cellbarcode_file"].notna(), "sample"].tolist()
        )