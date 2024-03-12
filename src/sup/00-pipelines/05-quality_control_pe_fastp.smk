rule qc_merged_pe:
    input:
        r1=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_1.fq.gz"),
        r2=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_2.fq.gz")
    output:
        out=os.path.join(config['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.trim.fq.gz"),
        html=os.path.join(config['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.html"),
        json=os.path.join(config['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.json")
    singularity:
        config['IMAGES']['fastp']
    params:
    threads: 8
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --merged_out {output.out} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            --cut_mean_quality 15 \
            --length_required 100 \
            --cut_front \
            --cut_tail \
            -q 15 \
            --merge \
            --include_unmerged \
            --detect_adapter_for_pe
        """