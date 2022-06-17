#!/usr/bin/env nextflow

/*

Pipeline for 16S and ITS data
to study co-abundance between bacteria fungi

*/

Channel
	.fromPath(params.runs_csv)
	.splitCsv(header:true)
	.branch{
		sra: it.run_id != "NA"
		local: it.run_id == "NA"
	}
	.set{source_runs}

process download_sra_runs {
	input:
	tuple val(sample_id), val(run_id), val(amplicon) from source_runs.sra.map{x -> [x.sample_id, x.run_id, x.amplicon]}

	tag "{\"sample_id\": \"$sample_id\", \"run_id\": \"$run_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("*.fq.gz")) into sra_raw_runs

	publishDir "${params.results_dir}/1-raw", mode: 'copy'
	conda "${params.conda_envs_dir}/bio.yml"
	maxForks 5

	script:
	"""
	grabseqs sra $run_id
	
	if [[ -e ${run_id}_1.fastq.gz ]] && [[ -e ${run_id}_2.fastq.gz ]] ; then
		echo run is paired!
		seq 1 2 | xargs -i  mv ${run_id}_{}.fastq.gz ${sample_id}_{}.${amplicon}.raw.fq.gz 

	elif [[ -e ${run_id}.fastq.gz ]] ; then
		echo run is unpaired!
		mv ${run_id}.fastq.gz ${sample_id}.${amplicon}.raw.fq.gz 
	fi
	"""
}

process copy_local_runs {
	input:
	tuple val(sample_id), val(library_layout), val(amplicon), val(sample_alias) from source_runs.local.map{x -> [x.sample_id, x.library_layout, x.amplicon, x.sample_alias]}

	tag "{\"sample_id\": \"$sample_id\", \"sample_alias\": \"$sample_alias\" , \"library_layout\": \"$library_layout\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("*.fq.gz")) into local_raw_runs

	publishDir "${params.results_dir}/1-raw", mode: 'copy'
	
	script:
	"""
	if [ "$library_layout" == "PAIRED" ]; then
		seq 1 2 | xargs -i cp -a ${params.raw_reads_dir}/${sample_alias}_{}.${amplicon}.raw.fq.gz \
			${sample_id}_{}.${amplicon}.raw.fq.gz
	else
		 cp -a ${params.raw_reads_dir}/${sample_alias}.${amplicon}.raw.fq.gz \
			${sample_id}.${amplicon}.raw.fq.gz
	fi
	"""
}

local_raw_runs
	.mix(sra_raw_runs)
	.into{raw_runs; fastqc_raw_runs}

raw_runs
	.branch{
		paired: it[2].size() == 2
			return [it[0], it[1], it[2][0], it[2][1]] // unnest pairs
		single: true
	}
	.set{library_layout_raw_runs}

process merge_paired_reads {
	input:
	tuple val(sample_id), val(amplicon), file("1.fq.gz"), file("2.fq.gz") from library_layout_raw_runs.paired

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("${sample_id}.${amplicon}.merged.fq.gz")) into paired_merged_runs

	publishDir "${params.results_dir}/2-merge", mode: 'copy'
	conda "${params.conda_envs_dir}/bio.yml"

	script:
	"""
	NGmerge -1 1.fq.gz -2 2.fq.gz -o ${sample_id}.${amplicon}.merged.fq.gz
	"""
}

process merge_single_reads {
	input:
	tuple val(sample_id), val(amplicon), file("raw.fq.gz") from library_layout_raw_runs.single

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("${sample_id}.${amplicon}.merged.fq.gz")) into single_merged_runs

	publishDir "${params.results_dir}/2-merge", mode: 'copy'

	script:
	"""
	ln raw.fq.gz ${sample_id}.${amplicon}.merged.fq.gz
	"""
}

paired_merged_runs
	.mix(single_merged_runs)
	.into{merged_runs; fastqc_merged_runs}

process run_after_merge_qc {
	input:
	tuple val(sample_id), val(amplicon), file("merged.fq.gz") from merged_runs

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("${sample_id}.${amplicon}.qc.fq.gz")) into run_after_merge_qc_output

	publishDir "${params.results_dir}/3-qc", mode: 'copy'
	conda "${params.conda_envs_dir}/bio.yml"

	script:
	"""
	# build default command
	cat > cmd.sh <<- EOF
	trimmomatic SE \
		merged.fq.gz \
		${sample_id}.${amplicon}.qc.fq.gz \
		ILLUMINACLIP:${params.adapters_path}/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
		LEADING:20 TRAILING:20 MINLEN:100 \
		-trimlog ${sample_id}.${amplicon}.qc.log
	EOF

	# try default command
	bash cmd.sh ||
		# try given phred score format
		(
			echo " -phred33" >> cmd.sh
			# remove line breaks
			cat cmd.sh | tr -d '\n' > cmd2.sh
			mv cmd2.sh cmd.sh

			bash cmd.sh
		) ||
		# fix file
		(
			echo $PATH > path

			fix_fastq.py -i merged.fq.gz -o fixed.fq.gz
			sed -i 's/merged.fq.gz/fixed.fq.gz/g' cmd.sh
			bash cmd.sh
		)
	"""
}

run_after_merge_qc_output
	.into{qc_runs; fastqc_qc_runs}

fastqc_raw_runs
	.mix(fastqc_merged_runs)
	.mix(fastqc_qc_runs)
	.map{x -> x[2]} // pluck file(s)
	.flatten()
	.set{fastqc_runs}

process run_fastqc {
	input:
	file fastq_file from fastqc_runs

	tag "{\"fastq_file\": \"$fastq_file\"}"	

	output:
	file "*_fastqc.{zip,html}" into fastqc_output

	conda "${params.conda_envs_dir}/bio.yml"

	script:
	"""
	fastqc -q $fastq_file
	"""
}


fastqc_output
	.flatten()
	// file path to directory
	.map{ x -> x.toString() - ~/[^\/]+$/}
	.unique()
	.ifEmpty([])
	.set{multiqc_input}

process run_multiqc {
	input:
	val fastqc_dir from multiqc_input

	tag "{\"fastq_dir\":\"$fastqc_dir\"}"

	output:
	file "multiqc_data/multiqc_fastqc.txt" into multiqc_output

	conda "${params.conda_envs_dir}/bio.yml"

	script:
	"""
	multiqc -m fastqc ${fastqc_dir}
	"""
}

process dump_multiqc_fastqc {
	//BUG: multiqc does not work on all directories at once
	//WORKARROUND: do it on per sample basis and only publish fastqc txts per sample
	input:
	file "run_multiqc_fastqc.txt" from multiqc_output.collectFile(newLine: true)

	output:
	file "multiqc_fastqc.txt"

	publishDir "${params.results_dir}/multiqc", mode: 'copy'

	script:
	"""
	cat run_multiqc_fastqc.txt >> multiqc_fastqc.txt
	"""
}

process import_qiime2 {
	input:
	tuple(val(sample_id), val(amplicon), val(qc_fq_file)) from qc_runs

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("qc.seqs.qza")) into qiime2_imported_runs

	conda "${params.conda_envs_dir}/qiime2-2020.2-py36-linux-conda.yml"

	script:
	"""
	echo -e "sampleid\tabsolute-filepath\n${sample_id}.${amplicon}\t${qc_fq_file}" > manifest.tsv

	qiime tools import \
		--type 'SampleData[SequencesWithQuality]' \
		--input-path manifest.tsv \
		--output-path qc.seqs.qza \
		--input-format SingleEndFastqManifestPhred33V2
	"""
}

process dereplicate_qiime2 {
	input:
	tuple(val(sample_id), val(amplicon), val(qc_qza_file)) from qiime2_imported_runs

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	

	output:
	tuple(val(sample_id), val(amplicon), file("dereplicated.seqs.qza"),  file("dereplicated.table.qza")) into qiime2_dereplicated_runs

	conda "${params.conda_envs_dir}/qiime2-2020.2-py36-linux-conda.yml"

	script:
	"""
	qiime vsearch dereplicate-sequences \
		--i-sequences $qc_qza_file \
		--o-dereplicated-table dereplicated.table.qza \
		--o-dereplicated-sequences dereplicated.seqs.qza 
	"""
}

qiime2_dereplicated_runs
	.branch {
		fun_its: it[1] == "fun_its"
		bac_16s: it[1] == "bac_16s"
	}
	.set{amplicon_qiime2_dereplicated_runs}

process pick_otus_bac_16s_qiime2 {
	input:
	tuple(val(sample_id), val(amplicon), file("dereplicated.seqs.qza"), file("dereplicated.table.qza")) from amplicon_qiime2_dereplicated_runs.bac_16s
	val perc_identity from params.bac_16s_perc_identity
	val db_taxonomy_qza from params.bac_16s_taxonomy_qza
	val db_seqs_qza from params.bac_16s_seqs_qza
	
	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"	
	
	output:
	tuple(val(sample_id), val(amplicon), file("clustered.table.qza"), file("clustered.seqs.qza"), file("clustered.unmatched.qza")) into bac_16s_qiime2_picked_otus_runs

	conda "${params.conda_envs_dir}/qiime2-2020.2-py36-linux-conda.yml"

	script:
	"""
	qiime vsearch cluster-features-closed-reference \
		--i-table dereplicated.table.qza \
		--i-sequences dereplicated.seqs.qza \
		--i-reference-sequences $db_seqs_qza \
		--p-perc-identity $perc_identity \
		--o-clustered-table clustered.table.qza \
		--o-clustered-sequences clustered.seqs.qza \
		--o-unmatched-sequences clustered.unmatched.qza \
		--p-threads ${task.cpus}
	"""
}

process pick_otus_fun_its_qiime2 {
	input:
	tuple(val(sample_id), val(amplicon), file("dereplicated.seqs.qza"), file("dereplicated.table.qza")) from amplicon_qiime2_dereplicated_runs.fun_its
	val perc_identity from params.fun_its_perc_identity
	val db_taxonomy_qza from params.fun_its_taxonomy_qza
	val db_seqs_qza from params.fun_its_seqs_qza

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"

	output:
	tuple(val(sample_id), val(amplicon), file("clustered.table.qza"), file("clustered.seqs.qza"), file("clustered.unmatched.qza")) into fun_its_qiime2_picked_otus_runs

	conda "${params.conda_envs_dir}/qiime2-2020.2-py36-linux-conda.yml"

	script:
	"""
	qiime vsearch cluster-features-closed-reference \
		--i-table dereplicated.table.qza \
		--i-sequences dereplicated.seqs.qza \
		--i-reference-sequences $db_seqs_qza \
		--p-perc-identity $perc_identity \
		--o-clustered-table clustered.table.qza \
		--o-clustered-sequences clustered.seqs.qza \
		--o-unmatched-sequences clustered.unmatched.qza \
		--p-threads ${task.cpus}
	"""
}

bac_16s_qiime2_picked_otus_runs
	.mix(fun_its_qiime2_picked_otus_runs)
	.set{qiime2_picked_otus_runs}

process export_otu_tables_qiime2 {
	input:
	tuple(val(sample_id), val(amplicon), file("clustered.table.qza"), file("clustered.seqs.qza"), file("clustered.unmatched.qza")) from qiime2_picked_otus_runs

	tag "{\"sample_id\": \"$sample_id\", \"amplicon\": \"$amplicon\"}"

	output:
	file "${sample_id}.${amplicon}.features.biom"

	conda "${params.conda_envs_dir}/qiime2-2020.2-py36-linux-conda.yml"
	publishDir "${params.results_dir}/4-features", mode: 'copy'

	script:
	"""
	qiime tools export \
		--input-path clustered.table.qza \
		--output-path .

	mv feature-table.biom  ${sample_id}.${amplicon}.features.biom
	"""
}

