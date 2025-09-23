from pathlib import Path


def get_cellbarcode_file(wildcards):
    return samples_df.loc[
        samples_df["sample"] == wildcards.sample, "cellbarcode_file"
    ].values[0]


def get_fastqs(wildcards):
    sample = wildcards.sample
    files = list(Path("trimmed").glob(f"{sample}*.fastq.gz"))
    print(files)
    return files


rule cellbarcode_file:
    input:
        get_cellbarcode_file,
    output:
        temp(".{sample}.cb_alias.txt"),
    shell:
        "cut -f1 {input} > {output}"


rule star_index:
    input:
        fasta=config["genome_dir"] + "/{genome}/{genome}.fa",
        gtf=config["genome_dir"] + "/{genome}/{genome}.annotation.gtf",
    output:
        directory=directory(config["genome_dir"] + "/{genome}/index/star/"),
        genome_parameters=config["genome_dir"]
        + "/{genome}/index/star/genomeParameters.txt",
    threads: 20
    log:
        "logs/star_index/{genome}.log",
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.directory} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 \
            &> {log}
        """


rule star_paired:
    input:
        r1="trimmed/{sample}_R1.fastq.gz",
        r2="trimmed/{sample}_R2.fastq.gz",
        index=lambda wildcards: expand(
            rules.star_index.output.directory,
            genome=samples_df.loc[
                samples_df["sample"] == wildcards.sample, "genome"
            ].values[0],
        ),
    output:
        bam="STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        log="STAR/{sample}/{sample}.Log.final.out",
    threads: 20
    params:
        prefix="STAR/{sample}/{sample}.",
    log:
        "logs/star/{sample}.log",
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.prefix} &> {log}
        """


rule star_solo:
    input:
        r1="trimmed/{sample}_R1.fastq.gz",
        r2="trimmed/{sample}_R2.fastq.gz",
        cell_barcodes=rules.cellbarcode_file.output,
        index=lambda wildcards: expand(
            rules.star_index.output.directory,
            genome=samples_df.loc[
                samples_df["sample"] == wildcards.sample, "genome"
            ].values[0],
        ),
    output:
        matrix="STAR/{sample}/{sample}_Solo.out/Gene/raw/matrix.mtx",
        features="STAR/{sample}/{sample}_Solo.out/Gene/raw/features.tsv",
        barcodes="STAR/{sample}/{sample}_Solo.out/Gene/raw/barcodes.tsv",
        log="STAR/{sample}/{sample}_Log.final.out",
    params:
        prefix="STAR/{sample}/{sample}_",
    log:
        "logs/star_solo/{sample}.log",
    threads: 20
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir  {input.index} \
            --readFilesIn {input.r2} {input.r1} \
            --readFilesCommand zcat \
            \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingThreadN {threads} \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outSAMmapqUnique 60 \
            --outSAMunmapped Within \
            --outFileNamePrefix {params.prefix} \
            \
            --outFilterMultimapNmax 1  \
            \
            --soloType CB_UMI_Simple \
            --soloCBstart 1 --soloCBlen 14 \
            --soloUMIstart 15 --soloUMIlen 14 \
            --clipAdapterType CellRanger4 \
            --soloUMIdedup 1MM_Directional \
            --soloCBmatchWLtype 1MM \
            --soloCellFilter None \
            --soloCBwhitelist {input.cell_barcodes} \
            --soloBarcodeReadLength 0 \
            --soloStrand Forward \
            --soloFeatures Gene \
            --quantMode GeneCounts &>> {log}
        """


rule matrix_to_tsv:
    input:
        matrix=rules.star_solo.output.matrix,
        features=rules.star_solo.output.features,
        barcodes=get_cellbarcode_file,
    output:
        tsv="counts/{sample}.tsv",
    wildcard_constraints:
        sample=r"(?!\.summaries/).*",
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
            --out {output.tsv} \
            &> {log}
        """


rule cellbarcode_summary:
    input:
        rules.matrix_to_tsv.output.tsv,
    output:
        "counts/.summaries/{sample}.cb_summary_mqc.tsv",
    log:
        "logs/cell_barcode_summary/{sample}.log",
    threads: 1
    shell:
        """
        set -euo pipefail
        python {workflow.basedir}/scripts/cellbarcode_summary.py \
            --count-table {input} \
            --output {output} \
            &> {log}
        """
