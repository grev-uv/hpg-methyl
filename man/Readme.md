
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
mapped read. This data is stored as individual CSV files, containing
the methylation status for each possible context (CpG, CHH, CHG and MUT).

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
