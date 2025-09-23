ruleorder: fastp_paired > fastp_single


rule sample_renaming:
    output:
        "trimmed/reports/.sample_renaming.tsv",
    run:
        with open(output[0], "w") as out_file:
            for _, row in samples_df.iterrows():
                sample = row["sample"]
                alias = row.get("alias") or sample
                out_file.write(f"{sample}\t{alias}\n")
                out_file.write(f"{sample}_R1\t{alias}\n")


rule fastp_paired:
    input:
        fq1="raw/{sample}_R1.fastq.gz",
        fq2="raw/{sample}_R2.fastq.gz",
    output:
        r1="trimmed/{sample}_R1.fastq.gz",
        r2="trimmed/{sample}_R2.fastq.gz",
        html="trimmed/reports/{sample}.html",
        json="trimmed/reports/{sample}.json",
    threads: 4
    log:
        "logs/fastp/{sample}.log",
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            -w {threads} &> {log}
        """


rule fastp_single:
    input:
        fq="raw/{sample}.fastq.gz",
    output:
        fq="trimmed/{sample}.fastq.gz",
        html="trimmed/reports/{sample}.html",
        json="trimmed/reports/{sample}.json",
    threads: 4
    log:
        "logs/fastp/{sample}.log",
    wildcard_constraints:
        sample=r"(?!.*_R[12]$)[A-Za-z0-9_.-]+",
    shell:
        """
        fastp -i {input.fq} -o {output.fq} \
        -h {output.html} -j {output.json} \
        -w {threads} &> {log}
        """


rule multiqc:
    input:
        fastp=expand(
            "trimmed/reports/{sample}.json", sample=samples_df["sample"].tolist()
        ),
        star=expand("STAR/{sample}/{sample}.Log.final.out", sample=bulk_samples),
        star_solo=expand("STAR/{sample}/{sample}_Log.final.out", sample=brbseq_samples),
        cellbarcode_summary=expand(
            "counts/.summaries/{sample}.cb_summary_mqc.tsv", sample=brbseq_samples
        ),
        sample_renaming="trimmed/reports/.sample_renaming.tsv",
    output:
        "multiqc_report.html",
    threads: 4
    log:
        "logs/multiqc/multiqc.log",
    shell:
        """
        set -euo pipefail
        {{
            rm -rv trimmed/reports/multiqc* || true
            multiqc {input} --config {workflow.basedir}/resources/multiqc_config.yaml --replace-names {input.sample_renaming} -o trimmed/reports/
            mv -v trimmed/reports/multiqc_report.html {output}
        }} &> {log}
        """
