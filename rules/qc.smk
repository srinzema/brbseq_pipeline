rule fastp:
    input: "raw/{sample}.done"
    output:    
        html="trimmed/reports/{sample}.html",
        json="trimmed/reports/{sample}.json",
    threads: 4
    log: "logs/fastp/{sample}.log"
    shell:
        """
        set -euo pipefail
        {{
            echo "Processing sample {wildcards.sample}"
            ls raw/{wildcards.sample}*

            FQ1="raw/{wildcards.sample}_R1.fastq.gz"
            FQ2="raw/{wildcards.sample}_R2.fastq.gz"
            SINGLE="raw/{wildcards.sample}.fastq.gz"

            REPORT_HTML=trimmed/reports/{wildcards.sample}.html
            REPORT_JSON=trimmed/reports/{wildcards.sample}.json

            if [[ -f $FQ1 && -f $FQ2 ]]; then
                fastp \
                    -i $FQ1 -I $FQ2 \
                    -o trimmed/{wildcards.sample}_R1.fastq.gz \
                    -O trimmed/{wildcards.sample}_R2.fastq.gz \
                    -h $REPORT_HTML  -j $REPORT_JSON \
                    -w {threads}
            elif [[ -f $SINGLE ]]; then
                fastp \
                    -i $SINGLE \
                    -o trimmed/{wildcards.sample}.fastq.gz \
                    -h $REPORT_HTML  -j $REPORT_JSON \
                    -w {threads}
            else
                echo "Error: no input FASTQs found for {wildcards.sample}"
                exit 1
            fi
        
            touch {output}
        }} &> {log}
        """

rule multiqc:
    input:
        expand("trimmed/reports/{sample}.html", sample=samples_df["alias"].tolist())
    output:
        "multiqc_report.html"
    threads: 4
    log: "logs/multiqc/multiqc.log"
    shell:
        """
        set -euo pipefail
        {{
            rm -rv trimmed/reports/multiqc* || true
            multiqc trimmed/reports/ -o trimmed/reports/
            mv -v trimmed/reports/multiqc_report.html {output}
        }} &> {log}
        """