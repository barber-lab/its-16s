rule get_raw_fastq_se:
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fastq")
    conda:
        os.path.join(config['workdir'], config['ENVIRONMENTS']['sra_tools'])
    threads: 4
    params:
         os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'])
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            -f \
            {wildcards.ID}
        """

rule fastqc_raw_se:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fastq"),
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_se'], "{ID}_fastqc.html")
    params:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_se'])
    singularity:
        config['IMAGES']['fastqc']
    threads: 4
    shell:
        """
        fastqc -t {threads} --outdir {params} {input}
        """

rule compress_fastq_se:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fastq")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fq.gz")
    threads: 4
    shell:
        """
        pigz -cp {threads} {input} > {output}
        """