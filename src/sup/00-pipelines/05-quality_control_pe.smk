rule merge_pe:
    input:
        r1=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_1.fq.gz"),
        r2=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_2.fq.gz"),
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fq.gz")
    threads: 4
    singularity:
        config['IMAGES']['ngmerge']
    shell:
        """
        NGmerge \
            -1 {input.r1} \
            -2 {input.r2} \
            -n {threads} \
            -o {output} \
            -u 41 \
            -g
        """

# rule fix_merged_fastq:
#     input:
#         os.path.join(config['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fq.gz")
#     output:
#         os.path.join(config['PROJ'], config['FOLDERS']['merged_pe_fixed'], "{ID}.fq.gz")
#     params: 
#         config['SCRIPTS']['fix_fq']
#     threads: 4
#     conda:
#         os.path.join(config['workdir'], "src/conda/python3.8.yml")
#     shell:
#         """
#         python {params} -i {input} -o {output}
#         """

rule Trimmomatic_PE:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fq.gz")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.fq.gz")
    params:
        "/usr/local/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
    singularity:
        config['IMAGES']['trimmomatic']
    threads: 4
    shell:
        """
        set +e
        trimmomatic SE \
            -threads {threads} \
            {input} \
            {output} \
            ILLUMINACLIP:{params}:2:30:10:2:keepBothReads \
            LEADING:20 \
            TRAILING:20 \
            MINLEN:100
        exitcode=$?

        if [ $exitcode -gt 0 ]; then
            # Retry with the additional argument
            trimmomatic SE \
                -phred33 \
                -threads {threads} \
                {input} \
                {output} \
                ILLUMINACLIP:{params}:2:30:10:2:keepBothReads \
                LEADING:20 \
                TRAILING:20 \
                MINLEN:100
            exitcode=$?
        fi

        set -e  # Reset the error handling
        exit $exitcode
        """

rule fastqc_trim_merged:
    input:
        os.path.join(config['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.fq.gz")
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_merg'], "{ID}_fastqc.html")
    params:
        os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_merg'])
    singularity:
        config['IMAGES']['fastqc']
    threads: 4
    shell:
        """
        fastqc -t {threads} --outdir {params} {input}
        """

# rule merge_pe:
#     input:
#         r1=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_1.fq.gz"),
#         r2=os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_2.fq.gz")
#     output:
#         os.path.join(config['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fq.gz")
#     threads: 4
#     singularity:
#         config['IMAGES']['bbtools']
#     shell:
#         """
#         bbmerge.sh \
#             in1={input.r1} \
#             in2={input.r2} \
#             out={output}
#         """