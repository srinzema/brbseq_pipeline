
def cellbarcode_file(wildcards):
    return samples_df.loc[samples_df["sample"] == wildcards.sample, "cellbarcode_file"].values[0]


rule star_index:
    input:
        fasta=config["genome_dir"] + "/{genome}/{genome}.fa",
        gtf=config["genome_dir"] + "/{genome}/{genome}.annotation.gtf"
    output: 
        directory(config["genome_dir"] + "/{genome}/index/star/")
    threads: 20
    log:
        "logs/star_index/{genome}.log"
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 \
            &> {log}
        """


rule star_solo:
    input: 
        r1="trimmed/{sample}_R1.fastq.gz",
        r2="trimmed/{sample}_R2.fastq.gz",
        cell_barcodes=cellbarcode_file,
        index=lambda wildcards: expand(
            rules.star_index.output, 
            genome=samples_df.loc[samples_df["sample"] == wildcards.sample, "genome"].values[0]
        )
    output:
        matrix="STAR/{sample}/{sample}_Solo.out/Gene/raw/matrix.mtx",
        features="STAR/{sample}/{sample}_Solo.out/Gene/raw/features.tsv",
        barcodes="STAR/{sample}/{sample}_Solo.out/Gene/raw/barcodes.tsv",
        log="STAR/{sample}/{sample}_Log.final.out",
    params:
        prefix="STAR/{sample}/{sample}_"
    log:
        "logs/star_solo/{sample}.log"
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
        barcodes=rules.star_solo.output.barcodes,
    output:
        tsv="counts/{sample}.tsv"
    log:
        "logs/solo_to_tsv/{sample}.log"
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
