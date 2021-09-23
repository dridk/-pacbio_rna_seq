from glob import glob
from itertools import product 
import re
import hashlib

RAW_DATA ="run1.fastq"

HG19 = config["REFERENCE"]

PLOT_SCRIPT = workflow.basedir +  "/plot_transcript_count.py"
HASH_BED_SCRIPT = workflow.basedir +  "/hash_bed.py"
PLOT_HASH_BED = workflow.basedir +  "/plot_hash_bed.py"


def fasta_name(filename):
	with open(filename) as file:
		return [line[1:].strip() for line in file if line.startswith(">")]


def all_output():

	barcode = fasta_name(config["BARCODE"])
	gene = [re.sub(r"_.", "", i) for i in fasta_name(config["PRIMERS"])]

	for i in product(gene,barcode):
		gene, barcode = i
		yield f"{gene}.{barcode}.hash.png"


def all_debarcoding():

	barcode = fasta_name(config["BARCODE"])
	for i in barcode:
		yield f"debarcoding.{i}--{i}.fastq"

print(all_output())

rule everything:
	input:
		list(all_output()),



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


rule debarcoding:
	input:
		RAW_DATA
	output:
		list(all_debarcoding())
	shell:
		"lima -W 300 {input} {config[BARCODE]} debarcoding.fastq --ccs --split-named --single-side --dump-clips"

# rule extract_barcode_and_trim:
# 	input:
# 		RAW_DATA
# 	output:
# 		"raw.{barcode}.fastq"
# 	shell:
# 		"""
# 		BARCODE=$(seqkit grep {config[BARCODE]} -p {wildcards.barcode} |seqkit seq -s)
# 		seqkit grep -s -p $BARCODE {input} > {output}
# 		"""


rule extract_amplicon:
	input: 
		"debarcoding.{barcode}--{barcode}.fastq"
	output:
		"{gene}.{barcode}.fastq"
	log:
		"{gene}.{barcode}.amplicon.stats"
	shell:
		"""
		FORWARD=$(seqkit grep {config[PRIMERS]} -p {wildcards.gene}_F|seqkit seq -s)
		REVERSE=$(seqkit grep {config[PRIMERS]} -p {wildcards.gene}_R|seqkit seq -s)
		seqkit amplicon {input} -F $FORWARD -R $REVERSE > {output}
		"""



rule read_distribution:
	input:
		"{gene}.{barcode}.fastq"
	output:
		"{gene}.{barcode}.reads.txt"
	shell:
		"""
		seqkit seq -s {input} | awk '{{print length($0)}}' > {output}
		"""

rule fastqc:
	input:
		"{gene}.{barcode}.fastq"
	output:
		"{gene}.{barcode}_fastqc.html"
	shell:
		"fastqc {input}"

# LOOK FOR ANTI SENS 
# MOTIF=$(seqkit grep primers.fa -p {wildcards.gene}|seqkit seq -srp)
# seqkit grep {input} -s -p $MOTIF | seqkit seq -rp > {output}

		# MOTIF=$(seqkit grep primers.fa -p {wildcards.gene}|seqkit seq -s)
		# seqkit grep {input} -s -p $MOTIF > {output}


rule minimap2:
	input:
		"{gene}.{barcode}.fastq"
	output:
		"{gene}.{barcode}.sam"
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
		"{gene}.{barcode}.bam"
	output:
		"{gene}.{barcode}.bed"
	shell:
		"bamToBed -bed12 -i {input} > {output}"

rule bed2hash:
	input:
		"{gene}.{barcode}.bed"
	output:
		"{gene}.{barcode}.hash.bed"
	shell:
		"python {HASH_BED_SCRIPT} {input} {output}"

rule plot_hashbed:
	input:
		"{gene}.{barcode}.hash.bed"
	output:
		"{gene}.{barcode}.hash.png"
	shell:
		"python {PLOT_HASH_BED} {input} {output}"


rule unique_hash:
	input:
		"{gene}.{barcode}.hash.bed"
	output:
		"{gene}.{barcode}.hash.unique"
	shell:
		"cat {input}|cut -f8|sort|uniq -c|sort -k1 -nr > {output}"

# rule uniqbed:
# 	input:
# 		"{gene}.{barcode}.bed"
# 	output:
# 		"{gene}.transcript.{barcode}.unique.bed"
# 	shell:
# 		"cat {input}|awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,\"name\",$5,$6,$7,$8,$9,$10,$11,$12}}'|sort |uniq -c|sort -k1 -nr|awk 'BEGIN{{OFS=\"\t\"}}{{print $2,$3,$4,$1, $6,$7,$8,$9,$10,$11,$12,$13}}' > {output}"


rule plotuniqbed:
	input:
		"{gene}.transcript.{barcode}.unique.bed"
	output:
		"{gene}.transcript.{barcode}.dist.png"
	shell:
		"python {PLOT_SCRIPT} {input} 10"
	

rule uniqbed50:
	input:
		"{gene}.transcript.{barcode}.unique.bed"
	output:
		"{gene}.transcript.{barcode}.min50.bed"
	shell:
		"cat {input}|awk '$4 > 50 {{ print $0}}' > {output}"
