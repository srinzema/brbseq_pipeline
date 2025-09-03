rule download_sra:
    output: temp("raw/{sample}.done")
    params: 
        accession=lambda wc: samples_df.loc[samples_df["alias"] == wc.sample, "sample"].values[0]
    threads: 4
    log:
        "logs/download/{sample}.log"
    shell:
        """
        set -euo pipefail
        {{
            fasterq-dump --split-files --threads {threads} {params.accession} -O raw/
            
            # paired-end if both _1 and _2 exist
            if [[ -f raw/{params.accession}_1.fastq && -f raw/{params.accession}_2.fastq ]]; then
                pigz -c raw/{params.accession}_1.fastq > raw/{wildcards.sample}_R1.fastq.gz
                pigz -c raw/{params.accession}_2.fastq > raw/{wildcards.sample}_R2.fastq.gz
            # otherwise single-end
            elif [[ -f raw/{params.accession}_1.fastq ]]; then
                pigz -c raw/{params.accession}_1.fastq > raw/{wildcards.sample}.fastq.gz
            else
                echo "No FASTQ found for {wildcards.sample}" >&2
                exit 1
            fi
            
            rm -fv raw/{params.accession}*.fastq

            touch {output}
        }} &> {log}
        """