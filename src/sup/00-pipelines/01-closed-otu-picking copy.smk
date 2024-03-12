workdir: config['workdir']
import pandas as pd
import os

# Functions to generate file names
def trim_merged_files(wildcards):
    pe_tbl = pd.read_csv(os.path.join(config['FOLDERS']['PROJ'], config['REFERENCE']['fq_pe']), header = 0, names = ["Run"], sep = "\t")
    pe_ids = list(pe_tbl["Run"])
    global PE
    PE = expand(os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.fastq"), ID = pe_ids)
    return PE

rule targets:
    input:
        trim_merged_files,
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['samples_artifact'], "samples.qza"),
        expand(os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['abundance'], "feature-table-abundances.tsv"))

rule merge_pe:
    input:
        r1=os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_1.fastq"),
        r2=os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_2.fastq")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fastq")
    threads: 12
    singularity:
        "https://depot.galaxyproject.org/singularity/ngmerge:0.3--ha92aebf_1"
    shell:
        """
        NGmerge \
            -1 {input.r1} \
            -2 {input.r2} \
            -n {threads} \
            -o {output}
        """


# ADD PARAMETERS FROM THE PAPER
rule qc_merged_pe:
    input:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['merged_pe'], "{ID}.fastq")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['trim_merg'], "{ID}.fastq")
    params:
        "/usr/share/trimmomatic/TruSeq3-PE.fa"
    singularity:
        "docker://biocontainers/trimmomatic:v0.38dfsg-1-deb_cv1"
    threads: 12
    shell:
        """
        TrimmomaticSE \
            -threads {threads} \
            {input} \
            {output} \
            ILLUMINACLIP:{params}:2:30:10:2:keepBothReads \
            LEADING:20 \
            TRAILING:20 \
            MINLEN:100
        """

rule create_manifest:
    input:
        trim_merged_files
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['manifest'], "manifest.tsv")
    run:
        with open(output[0], "w") as manifest:
            manifest.write("sample-id\tabsolute-filepath\n")

            for file_path in input:
                # Extract the sample ID and file path from the input file path
                sample_id = file_path.split("/")[-1].split(".")[0]
                full_path = f"{config['workdir']}/{file_path}"
                manifest.write(f"{sample_id}\t{full_path}\n")

# Import samples
rule qiime_import_data:
    input:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['manifest'], "manifest.tsv")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['samples_artifact'], "samples.qza")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path {input} \
            --output-path {output} \
            --input-format SingleEndFastqManifestPhred33V2
        """

rule qiime_import_16S:
    input:
        os.path.join(config['REFERENCE']['silva_unite'], "silva_seq.fna")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['reference_artifact'], "silva_seq.qza")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input} \
            --output-path {output}
        """

rule qiime_import_ITS:
    input:
        os.path.join(config['REFERENCE']['silva_unite'], "unite_seq.fna")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['reference_artifact'], "unite_seq.qza")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input} \
            --output-path {output}
        """

rule qiime_dereplicate:
    input:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['samples_artifact'], "samples.qza")
    output:
        tbl = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['dereplicated'], "samples.tbl.qza"),
        seq = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['dereplicated'], "samples.seq.qza")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime vsearch dereplicate-sequences \
            --i-sequences {input} \
            --o-dereplicated-table {output.tbl} \
            --o-dereplicated-sequences {output.seq}
        """

rule qiime_closed_reference:
    input:
        ref = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['reference_artifact'], "{db}_seq.qza"),
        tbl = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['dereplicated'], "samples.tbl.qza"),
        seq = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['dereplicated'], "samples.seq.qza")
    output:
        tbl = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['otu_closed'], "samples_{db}.tbl.qza"),
        seq = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['otu_closed'], "samples_{db}.seq.qza"),
        unm = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['otu_closed'], "samples_{db}.unmatch.qza")
    threads: 48
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime vsearch cluster-features-closed-reference \
            --p-threads {threads} \
            --i-table {input.tbl} \
            --i-sequences {input.seq} \
            --i-reference-sequences {input.ref} \
            --p-perc-identity 0.97 \
            --o-clustered-table {output.tbl} \
            --o-clustered-sequences {output.seq} \
            --o-unmatched-sequences {output.unm}
        """

rule qiime_export:
    input:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['otu_closed'], "samples_{db}.tbl.qza")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['export'], "{db}-feature-table.biom")
    params:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['export'])
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {params}
        
        mv {params}/feature-table.biom {output}
        """

rule add_taxonomy2biom:
    input:
        biom = os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['export'], "{db}-feature-table.biom"),
        tax = os.path.join(config['REFERENCE']['silva_unite'], "{db}_tax.txt")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['add_tax'], "{db}-feature-table.biom")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        biom add-metadata \
            -i {input.biom} \
            -o {output} \
            --observation-metadata-fp {input.tax} \
            --observation-header OTUID,TAXONOMY \
            --sc-separated taxonomy
        """

rule abundance_from_biom:
    input:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['add_tax'], "{db}-feature-table.biom")
    output:
        os.path.join(config['FOLDERS']['PROJ'], config['FOLDERS']['abundance'], "{db}-feature-table-abundances.tsv")
    singularity:
        config['IMAGES']['qiime2']
    shell:
        """
        biom convert \
            -i {input} \
            -o {output} \
            --to-tsv
        """