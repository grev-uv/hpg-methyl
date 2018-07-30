
HPG-Methyl Manual
=================

A PDF version of the user manual can be downloaded from [here](manual.pdf).

# Running

Run HPG-Methyl with:

```hpg-methyl <mode> <options>```

Where mode can be either `build-index`, which must be called before doing any
mappings to create the BWT index, and `bs`, used to map sequences to the genome
and extract the methylation status.

## BS mode

The `bs` mode is used to map sequenced reads with the reference genome and,
optionally, extract the methylation context information from the alignments.
The mandatory command line options for the mode are:

* `-i` / `--bwt-index`: Path to the directory containing the BWT index created
  previously with the `build-index` option.
* `-f` / `--fq` / `--fastq`: Path to the FASTQ sequence file containing the
  sequenced reads.
* `--cpu-threads`: The number of CPU threads to use. HPG-Methyl can fully exploit
  all the processing cores available in your computer platform, improving the
  performance.

## Burrows-Wheeler Index generation

To be able to run HPG-Methyl, the index used to access the reference genome must
be created first. This is a preprocess step and must be done only once for each
reference genome. The process can be tuned with the following options:

* `-g` / `--ref-genome`: The path to the FASTA file containing the reference
  genome.
* `-i` / `--bwt-index`: The path to the directory where the BWT index will
be stored.
* `-r` / `--index-ratio`: Index compression ratio. Default value is 10.
* `--bs-index`: Enable this option to create a BWT index compatible with the
  methylation status extraction process.



## Methylation status

When running with the `--write-mcontext` command line option, the
application performs an analysis of the methylation status for every
mapped read.

This data is stored as individual CSV files, containing the methylation
status for each possible context (CpG, CHH, CHG and MUT), and in the
optional tags of the BAM file alignments, following the naming used by
[Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/).

### CSV File Format

The output CSV files have the following fields:

* **Query name:** The name of the read containing the methylated C as present
  on the FASTQ file.
* **Status:** **+** if methylated, **-** if unmethylated, **.** if irrelevant.
* **Chromosome:** Chromosome index of the read.
* **Start:** Starting position of the methylated C in the chromosome.
* **Context:** Methylation context of the methylated C, following the same
  naming convention as Bismark (z/Z for un/methylated CpG context, x/X for
  un/methylated CHG context, h/H for un/methylated CHH context and u/U for
  un/methylated MUT context).

Each field are divided in columns using tabulations, and in rows using lines.
An example entry of the methylation output is:


`1_89741628_1_0_1_0_0_0:1:0_0:0:0_1a2 -	0	89741698	h	0`

| Query Name | Status | Chromosome index | Start | Context | Strand |
|------------|--------|------------|-------|---------|--------|
|1_89741628_1_0_1_0_0_0:1:0_0:0:0_1a2|-|0|89741698|Unmeth CHH|Positive|

### BAM Optional Tag Format

The alignments in the output BAM file have the following optional tags:

| Tag name | Type | SAM Type | Description | Possible values |
|----------|------|----------|-------------|-----------------|
|   AS     |Number| i        | Alignment score. |0 or greater |
|   NH     |Number| i        | Number of reported alignments containing the query in the current record.| 0 or greater |
|   NM     |Number| i        | Edit distance to the reference, including ambiguous bases but excluding clipping.| 0 or greater |
|   XM     |String| Z        | Per-base methylation context | **z** / **Z**, **x** / **X**, **h** / **H**, **u** / **U**, **.** (dot)|
|   XG     |Number| i        | Alignment conversion state | **CT**, **GA** |
|   XR     |Number| i        | Read conversion state | **CT**, **GA** |
|   ZM     |Number| i        | Number of methylated C's in the alignment | 0 ~ Sequence length |

Tags can be extracted using `samtools` with the `view` command, following is
an example query of a methylated SAM alignment:


```
Example 256 1 234085256 234 2H34M2D39M * 0 0 GGGAATTAA(···)TTAGTG a!!!(···)!!!! AS:i:0 NH:i:5 NM:i:1 ZM:i:6 XG:Z:CT XR:Z:CT XM:Z:.......hh....(···)...M.......
```

## Mapped read filtering

HPG-Methyl will write to the BAM file all the correctly mapped alignments. If
only certain alignments are needed, the output can be filtered to only record
certain alignments in the BAM file using the following options:

* `--report-all`: Store all the aligned reads.
* `--report-best`: If there are mapped reads matching to the same location on the
  reference genome, only store those with the highest score.
* `--report-n-first`: If there are mapped reads matching to the same location on the
  reference genome, store the first `<n>` hits found, regardless of their quality.
* `--report-n-best`: If there are mapped reads matching to the same location on the
  reference genome, store the `<n>` best hits found for the location.


## Paired-end mode

HPG-Methyl is able to process paired-end sequenced reads.
In order to do that, the `--paired-mode` command line option must be set to 1 (by default is 0, single-end mode).
The paired-end sequences must be separated in two differents fastq files, but only one bam file with the mapped alignments will be created.

If paired-end mode is set, a mandatory command line `-j` / `--fq2` / `--fastq2` is required, in order to set the path to the second FASTQ sequence file.

Following is an example of use:

hpg-methyl bs -i index-directory --paired-mode 1 -f fastq-file-path-1 -j fastq-file-path-2 -o output-directory --cpu-threads thread-count






