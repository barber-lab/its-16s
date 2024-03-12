rule get_raw_fastq_pd:
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fastq")
    conda:
        os.path.join(config['workdir'], config['ENVIRONMENTS']['sra_tools'])
    threads: 4
    params:
         os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'])
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            -f \
            {wildcards.ID}
        """

rule fastqc_raw_pe:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fastq"),
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_pe'], "{ID}_{PE}_fastqc.html")
    params:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_pe'])
    singularity:
        config['IMAGES']['fastqc']
    threads: 4
    shell:
        """
        fastqc -t {threads} --outdir {params} {input}
        """

rule compress_fastq_pe:
    input:
        link = rules.fastqc_raw_pe.output,
        gz = os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fastq")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fq.gz")
    threads: 4
    shell:
        """
        pigz -cp {threads} {input.gz} > {output}
        """