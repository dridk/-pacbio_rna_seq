This pipeline was created as part of the [GOLD project](This pipeline was created as part of the GOld project).


## Installation
#### Dependencies 
* [python >= 3.9 ](https://www.python.org/downloads)
   - seaborn
   - pandas
   - matplotlib
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [lima](https://lima.how/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [minimap2](https://lh3.github.io/minimap2/)
* [samtools](http://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/) 

#### Install environment from conda 

```bash
conda env create -n gold -f env.yaml
````

## Usage 

#### Clone the repository

```bash
git clone git@github.com:dridk/pacbio_rna_seq.git
```

#### Edit config.yaml
- ```FASTQ``` The Fastq file path generated by PacBio Sequencing 
- ```BARCODE``` The Fasta file path describing barcodes used by lima for demultiplexing ( see example in repository ) 
- ```PRIMERS``` The Fasta file describing primers used for PacBio amplicon sequencing ( see example in repository ) 
- ```REFERENCE``` The fasta reference file used by minimap2 for alignement ( e.g: hg19.fa ) 


#### Run the pipeline 

Put ```your_file.fastq``` generated by PacBio in the same folder than *config.yaml* and run the following command. 
You can edit how many threads you want to use with ```--cores``` option.

```
snakemake -Fp --cores 10 --configfile config.yaml 
```

## Output 

The pipeline will generate one file per barcode and amplicon. 
For instance HBB.bc1022.bam contains aligned reads from HBB amplicon and bc1022 barcode identifer.

- debarcoding.{barcode}--{barcode}.fastq : Demultiplexed reads 
- ```{amplicon}.{barcode}.fastq```  : Transcripts reads
- ```{amplicon}.{barcode}.bam```  : Aligned transcripts Reads 
- ```{amplicon}.{barcode}.bed```  : Transcripts structures as a bed file 
- ```{amplicon}.{barcode}.hash.bed```  : Transcripts structures as a bed file with a unique ID to identify the transcript
- ```{amplicon}.{barcode}.hash.png```  : Distribution plot of transcripts
- ```cluster.{amplicon}.png```  : Transcripts abundance heatmap 











