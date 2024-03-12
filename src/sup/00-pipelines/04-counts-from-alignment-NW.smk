rule indexing_bam:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam.bai")
    singularity:
        config['IMAGES']['samtools']
    threads: 8
    shell:
        """
        samtools index {input} {output}
        """

rule report_alignment_summary:
    input:
        set_order=rules.indexing_bam.output,
        bam=os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['stats'], "{ID}_{BASE_NAME}.txt")
    singularity:
        config['IMAGES']['samtools']
    threads: 1
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """