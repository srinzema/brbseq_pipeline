def get_cellbarcode_file(wildcards):
    return samples_df.loc[
        samples_df["sample"] == wildcards.sample, "cellbarcode_file"
    ].values[0]


rule matrix_to_tsv:
    input:
        matrix=rules.star_solo.output.matrix,
        features=rules.star_solo.output.features,
        barcodes=rules.star_solo.output.barcodes,
        aliases=get_cellbarcode_file,
    output:
        tsv="counts/brb_{sample}.tsv",
    log:
        "logs/solo_to_tsv/{sample}.log",
    threads: 1
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        python {workflow.basedir}/scripts/matrix_to_tsv.py \
            --matrix {input.matrix} \
            --features {input.features} \
            --barcodes {input.barcodes} \
            --aliases {input.aliases} \
            --out {output.tsv} \
            &> {log}
        """

rule bam_to_tsv:
    input:
        bam=rules.star_paired.output.bam,
        gtf=lambda wildcards: expand(
            config["genome_dir"] + "/{genome}/{genome}.annotation.gtf",
            genome=samples_df.loc[
                samples_df["sample"] == wildcards.sample, "genome"
            ].values[0],
        )
    output:
        tsv="counts/bulk_{sample}.tsv",
    log:
        "logs/bam_to_tsv/{sample}.log",
    threads: 1
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        featureCounts -a {input.gtf} -o {output.tsv} -T {threads} {input.bam} &> {log}
        """