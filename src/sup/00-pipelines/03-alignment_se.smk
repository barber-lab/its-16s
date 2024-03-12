rule bwa_mem2_generate_indexes:
    input:
        os.path.join(config['REFERENCE'], "{BASE_NAME}.fna"),
    output:
        os.path.join(config['PROJ'], config['GENOME_INDEX'],  "{BASE_NAME}.done"),
    singularity:
        config['IMAGES']['bwaMem2']
    params:
        idx_dir=os.path.join(config['PROJ'], config['GENOME_INDEX']),
        idx_base = os.path.join(config['PROJ'], config['GENOME_INDEX'], "{BASE_NAME}")
    threads: 48
    shell:
        """
        mkdir -p {params.idx_dir} &&
        bwa-mem2 index -p {params.idx_base} {input}
        touch {output}
        """

rule bwa_mem2_alignment_single:
    input:
        ref=os.path.join(config['PROJ'], config['GENOME_INDEX'],  "{BASE_NAME}.done"),
        r=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fq.gz"),
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['alngmt'], "{ID}_{BASE_NAME}.sam")
    params:
        idx_base = os.path.join(config['PROJ'], config['GENOME_INDEX'], "{BASE_NAME}")
    singularity:
        config['IMAGES']['bwaMem2']
    resources:
        mem_mb=100000
    threads: 8
    shell:
        """
        bwa-mem2 mem -t {threads} {params} {input.r} > {output}
        """

rule sam_to_bam:
    input:
        sam=os.path.join(config['PROJ'], config['FOLDERS']['alngmt'], "{ID}_{BASE_NAME}.sam")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['bam'], "{ID}_{BASE_NAME}.bam")
    singularity:
        config['IMAGES']['samtools']
    threads: 8
    shell:
        """
        samtools view -@ {threads} -hbS {input.sam} > {output}
        """

rule sorting_bam:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['bam'], "{ID}_{BASE_NAME}.bam")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam")
    singularity:
        config['IMAGES']['samtools']
    threads: 8
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        """