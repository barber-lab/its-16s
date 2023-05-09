# Get only bacterial sequences
seqkit grep -Inrp Bacteria SILVA_138.1_SSURef_tax_silva.fasta > silva138.1ssuRef_bacteria.fasta

# RNA to DNA
sed 's/U/T/g' silva138.1ssuRef_bacteria.fasta > silva138.1ssuRef_bacteria_dna.fasta

# Extracted hypervariable regions from Silva bacterial DNA using primers from https://doi.org/10.1101/2021.09.03.455391
perl /home/qi47rin/apps/in_silico_pcr/in_silico_PCR.pl -s silva138.1ssuRef_bacteria_dna.fasta -p ../02-16s/16s_primers.tsv > out.fa

# Remove duplicated sequences
seqkit rmdup -s -o clean.fa.gz silva138.1ssuRef_bacteria_16sHypervariableRegions.fa

# Rename duplicated headers
seqkit rename file.fasta
