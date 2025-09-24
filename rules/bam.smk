def get_cellbarcode_file(wildcards):
    return samples_df.loc[
        samples_df["sample"] == wildcards.sample, "cellbarcode_file"
    ].values[0]


rule index_solo:
    input:
        bam="STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        bai="STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/index_solo/{sample}.log",
    shell:
        "samtools index {input.bam} &> {log}"


rule cellbarcode_summary:
    input:
        bam="STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        bai="STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
        barcodes=get_cellbarcode_file,
    output:
        "counts/.summaries/{sample}.cb_summary_mqc.tsv",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/cell_barcode_summary/{sample}.log",
    threads: 6
    shell:
        """
        set -euo pipefail
        python {workflow.basedir}/scripts/cellbarcode_summary.py \
            --bam {input.bam} --cellbarcodes {input.barcodes} \
            --threads {threads} --output {output} \
            &> {log}
        """
