from glob import glob
from itertools import product 
import re


RAW_DATA ="run1.fastq"

HG19 = config["REFERENCE"]

def fasta_name(filename):
	with open(filename) as file:
		return [line[1:].strip() for line in file if line.startswith(">")]


def all_output():

	barcode = fasta_name(config["BARCODE"])
	gene = [re.sub(r"_.", "", i) for i in fasta_name(config["PRIMERS"])]

	for i in product(gene,barcode):
		gene, barcode = i
		yield f"{gene}.transcript.{barcode}.dist.png"




rule everything:
	input:
		list(all_output()),

rule extract_amplicon:
	input: 
		RAW_DATA
	output:
		"{gene}.merge.fastq"
	log:
		"{gene}.amplicon.stats"
	shell:
		"""
		FORWARD=$(seqkit grep primers.fa -p {wildcards.gene}_F|seqkit seq -s)
		REVERSE=$(seqkit grep primers.fa -p {wildcards.gene}_R|seqkit seq -s)
		seqkit seq {input}|seqkit amplicon -r -1000:1000 -f -F $FORWARD -R $REVERSE > {output}
		"""


# rule extract_amplicon_rc:
# 	input: 
# 		RAW_DATA
# 	output:
# 		"{gene}.amplicon.rc.fastq"
# 	log:
# 		"{gene}.amplicon.rc.log"
# 	shell:
# 		"""
# 		FORWARD=$(seqkit grep primers.fa -p {wildcards.gene}_F|seqkit seq -s)
# 		REVERSE=$(seqkit grep primers.fa -p {wildcards.gene}_R|seqkit seq -s)
# 		seqkit seq -p {input} |seqkit amplicon -F $FORWARD -R $REVERSE   > {output}
# 		"""

# rule merge_amplicon:
# 	input:
# 		"{gene}.amplicon.rc.fastq",
# 		"{gene}.amplicon.fastq"
# 	output:
# 		"{gene}.merge.fastq"
# 	shell:
# 		"cat {input}> {output}"


rule extract_barcode_and_trim:
	input:
		"{gene}.merge.fastq"
	output:
		"{gene}.merge.{barcode}.fastq"
	shell:
		"""
		BARCODE=$(seqkit grep barcode.fa -p {wildcards.barcode} |seqkit seq -s)
		seqkit grep -s -p $BARCODE {input} | seqkit subseq -r 2:-17 > {output}
		"""

rule read_distribution:
	input:
		"{gene}.merge.{barcode}.fastq"
	output:
		"{gene}.merge.{barcode}.reads.txt"
	shell:
		"""
		seqkit seq -s {input} | awk '{{print length($0)}}' > {output}
		"""

rule fastqc:
	input:
		"{gene}.merge.{barcode}.fastq"
	output:
		"{gene}.merge.{barcode}_fastqc.html"
	shell:
		"fastqc {input}"

# LOOK FOR ANTI SENS 
# MOTIF=$(seqkit grep primers.fa -p {wildcards.gene}|seqkit seq -srp)
# seqkit grep {input} -s -p $MOTIF | seqkit seq -rp > {output}

		# MOTIF=$(seqkit grep primers.fa -p {wildcards.gene}|seqkit seq -s)
		# seqkit grep {input} -s -p $MOTIF > {output}


rule minimap2:
	input:
		"{gene}.merge.{barcode}.fastq"
	output:
		"{gene}.merge.{barcode}.sam"
	threads:
		5
	shell:
		"minimap2 -t {threads} -ax splice:hq -uf {HG19} {input} > {output}"

rule sam2bam:
	input:
		"{filename}.sam"
	output:
		"{filename}.bam"
	shell:
		"samtools sort -O BAM {input} > {output}; samtools index {output}"


rule bam2bed:
	input:
		"{gene}.merge.{barcode}.bam"
	output:
		"{gene}.merge.{barcode}.bed"
	shell:
		"bamToBed -bed12 -i {input} > {output}"


rule uniqbed:
	input:
		"{gene}.merge.{barcode}.bed"
	output:
		"{gene}.transcript.{barcode}.unique.bed"
	shell:
		"cat {input}|awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,\"name\",$5,$6,$7,$8,$9,$10,$11,$12}}'|sort |uniq -c|sort -k1 -nr|awk 'BEGIN{{OFS=\"\t\"}}{{print $2,$3,$4,$1, $6,$7,$8,$9,$10,$11,$12,$13}}' > {output}"


rule plotuniqbed:
	input:
		"{gene}.transcript.{barcode}.unique.bed"
	output:
		"{gene}.transcript.{barcode}.dist.png"
	shell:
		"python plot_transcript_count.py {input} 10"
	

rule uniqbed50:
	input:
		"{gene}.transcript.{barcode}.unique.bed"
	output:
		"{gene}.transcript.{barcode}.min50.bed"
	shell:
		"cat {input}|awk '$4 > 50 {{ print $0}}' > {output}"
