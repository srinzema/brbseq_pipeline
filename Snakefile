import os
import pandas as pd
from pathlib import Path


configfile: "config.yaml"


samples_df = pd.read_csv(config["samples"], sep="\t", comment="#")
working_dir: config["outdir"]


# Set brbseq samples
for d in config["cellbarcodes"]:
    sample = d["sample"]
    cellbarcode_file = Path(d["file"])

    assert (
        cellbarcode_file.exists()
    ), f"Cell barcode file {cellbarcode_file} does not exist"
    samples_df.loc[samples_df["sample"] == sample, "cellbarcode_file"] = (
        cellbarcode_file
    )

bulk_samples = samples_df.loc[samples_df["cellbarcode_file"].isna(), "sample"].tolist()
brbseq_samples = samples_df.loc[
    samples_df["cellbarcode_file"].notna(), "sample"
].tolist()


include: "rules/download.smk"
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/bam.smk"


rule all:
    input:
        rules.multiqc.output,
        expand(rules.star_paired.output.bam, sample=bulk_samples),
        expand(
            rules.cellbarcode_summary.output,
            sample=samples_df.loc[
                samples_df["cellbarcode_file"].notna(), "sample"
            ].tolist(),
        ),
