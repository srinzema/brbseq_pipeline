
rule download_sra:
    output: temp("raw/{sample}.done")
    threads: 4
    log:
        "logs/download/{sample}.log"
    shell:
        """
        set -euo pipefail
        {{
            if [[ -f raw/{wildcards.sample}_R1.fastq.gz || -f raw/{wildcards.sample}.fastq.gz ]]; then
                echo "FASTQ files for {wildcards.sample} already exist, skipping download."
                touch {output}
                exit 0
            fi

            fasterq-dump --split-files --threads {threads} {wildcards.sample} -O raw/
            
            # paired-end if both _1 and _2 exist
            if [[ -f raw/{wildcards.sample}_1.fastq && -f raw/{wildcards.sample}_2.fastq ]]; then
                pigz -c raw/{wildcards.sample}_1.fastq > raw/{wildcards.sample}_R1.fastq.gz
                pigz -c raw/{wildcards.sample}_2.fastq > raw/{wildcards.sample}_R2.fastq.gz
            
            # otherwise single-end
            elif [[ -f raw/{wildcards.sample}_1.fastq ]]; then
                pigz -c raw/{wildcards.sample}_1.fastq > raw/{wildcards.sample}.fastq.gz
            
            else
                echo "No FASTQ found for {wildcards.sample}" >&2
                exit 1
            fi
            
            rm -fv raw/{wildcards.sample}*.fastq

            touch {output}
        }} &> {log}
        """


rule download_genome:
    output:
        fasta=config["genome_dir"] + "/{genome}/{genome}.fa",
        gtf=config["genome_dir"] + "/{genome}/{genome}.annotation.gtf"
    log:
        "logs/download/{genome}.log"
    threads: 8
    params:
        genome_dir=config["genome_dir"]
    conda:
        "../envs/download.yaml"
    shell:
        """
        set -euo pipefail
        {{
            echo "Downloading {wildcards.genome}"
            
            genomepy install {wildcards.genome} \
                -g {params.genome_dir} \
                -l {wildcards.genome} \
                -t {threads} \
                -a

        }} &>> {log}
        """
