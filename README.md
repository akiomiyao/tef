# TEF
## Transposable Element Finder - Detection of active transposable elements from NGS data  
Transposable Element Finder (TEF) is a general‑purpose, organism‑agnostic program
for detecting active transposable elements (TEs) from next‑generation sequencing (NGS) data
by comparing two samples, or one sample and reference genome.  

TEF detects transposition events without requiring prior knowledge of TE families or repeat libraries.
Instead, it identifies transposed transposon ends and their insertion positions
on a reference genome using sequence features derived directly from short‑read data.

## Concept and Method
Active transposable elements generate new insertions in the host genome.
Because individual transposition events are independent, the insertion sites of active TEs
differ between samples.  

Short‑read NGS data therefore contain fusion fragments combining genomic sequences
and transposon ends at newly inserted loci.  

Most transposable elements create target site duplications (TSDs) at the insertion site.
TEF detects both 5′ and 3′ transposon ends and determines their genomic insertion positions
by exploiting TSD information present in short‑read sequencing data.  

This comparative strategy enables direct detection of ongoing or recent transposition activity
without reliance on species‑specific TE annotations.

## Detection Methods

TEF supports two detection methods.

### Junction method (default)

- Detects transposition sites using junction reads.
- Requires a reference genome.
- Outputs genomic insertion positions.
- Generates verification files based on head- and tail-side alignments.

### TSD method

- Detects transpositions based on target site duplications (TSDs).
- Can be used with or without a reference genome.
- When no reference genome is provided, insertion positions on the genome
  are not determined.

Result files for the TSD method are the same as those for the junction method,
except that file names begin with `tsd_method`.
### Verification of insertion events

TEF extracts candidate insertion events based on junction reads and generates genotypes for each event.
From these genotypes, TEF performs read alignments to the genomic sequences
corresponding to the predicted transposition sites.

Only insertion events supported by both head- and tail-side alignments of the
transposable element are written to `junction_method.verify.a.b`.

## Scope and Applicability
TEF is not species‑specific and can be applied to any organism
for which NGS data from individuals, accessions,
or samples are available.  

Typical applications include:  
- Detection of active transposable elements from paired NGS datasets  
- Genome‑wide identification of TE insertion polymorphisms  
- Comparative and population‑scale analysis of TE activity  
- Pangenome studies and structural variation analyses  
- Evolutionary and functional genomics across diverse species

## Example Application: Rice (Oryza sativa)
TEF has been extensively applied to rice (Oryza sativa),
where it was used to generate a genome‑wide map of transposable element insertion polymorphisms
across multiple accessions.
This rice dataset provides a well‑characterized reference example demonstrating
how TEF can be used in pangenome analysis, evolutionary genomics,
and studies of domestication and breeding.

## System Requirements
TEF runs on Unix‑like systems (Linux and FreeBSD).  

Memory: 8–16 GB RAM (depending on dataset size)  
Storage: 1 TB or more recommended (SSD strongly recommended)  

For large genomes (e.g. human WGS data), 4 TB or more is recommended  

## Resources
Project web page:  
https://akiomiyao.github.io/tef/

## Usage

For example,  
```
perl tef.pl a=ttm2,ref=IRGSP1.0  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0  
```

Options are separated by commas WITHOUT spaces.

With the new feature, you can now omit the b parameter. When b is omitted, TEF will use the reference genome specified by ref as the second set of data. This means TEF will compare the NGS reads specified by a directly with the reference genome to detect transpositions.

TEF is implemented with junction and TSD methods. The default is the junction method.

If the tsd_size option is specified, TEF will run with the TSD method.
For example,  
```
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5  
```
The th option is threshold for filtering same TSD sequence.  
Default = 0.2  
If too many noises are detected, increase th value, *e.g.* th=0.7.  
```
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,th=0.7  
```
The th option is TSD method specific.  

Because temporary directories in targets grow to a huge size, they will be deleted at the end of the analysis. If debug=yes is added to the options, temporary directories will not be deleted.
  
sort_tmp=directory_for_sort is the option for specifying a temporary directory for the sort command. If sort_tmp is specified to a fast disk, e.g., SSD, sorting will be accelerated. If your SSD is not large enough, set the TEF directory in HDD space and specify sort_tmp in SSD space.

For Linux, the max_process is limited to the number of CPU cores. For other OS, the default process number is 4. If you add the max_process option, e.g., max_process=8, TEF will use 8 cores.
```
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5,th=0.7,option=clear,sort_tmp=/mnt/ssd/tmp,max_process=8
```
Additionally, if an error occurs and the program stops, execute ulimit -n 4096 before running this program. This will increase the limit on the number of files that can be opened at once.

## Input Files

- `read/`  
  Input FASTQ files (uncompressed or compressed: `.gz`, `.bz2`, `.xz`).  
  Both paired-end and single-end reads are supported.
  
## Output Files (Junction Method)

Result files are output into a directory specified by `a` and `b`.

### Main files

- `junction_method.genotype.a.b`  
  Intermediate genotype file. This is the central file for the junction method,
  and most other result files are generated from this file.

- `junction_method.summary.a.b`  
  Summary of detected transposition events.

- `junction_method.tepos.a.b`  
  Estimated genomic positions of transposable element insertions.

- `junction_method.a.b.vcf`  
  Detected insertion events in VCF format.

- `junction_method.verify.a.b`  
  Verified insertion events supported by both head- and tail-side alignments
  of the transposable element.  
  This file is intended for manual or downstream verification of insertion events.

- `junction_method.log`  
  Log file for the junction method.

### Supporting files and directories

- `count.20/`  
  Count data used to detect junctions.
  If additional analysis is required, keep this directory.

- `reads_containing_junction.a.b`  
  Reads containing junction sequences used for verification.

- `alignment.a.b`  
  Read alignment files generated to create the verify file.



## Getting tef.pl
Clone the repository:
```
git clone https://github.com/akiomiyao/tef.git   
```
for updating tef.pl  
```
cd tef  
git pull  
```
  
## Analysis by tef.pl  
First, create directories specified by a and b. 

For example,  
```
cd tef  
mkdir IRGSP1.0  # If tef.pl runs without ref directory, tef.pl makes ref directory and exits.
cd IRGSP1.0  
wget  https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz  
cd ..  
mkdir ttm2  
mkdir ttm2/read  
cp somewhere/SRR556173_1.fastq somewhere/SRR556173_2.fastq ttm2/read  
mkdir ttm5  
mkdir ttm5/read  
cp somewhere/SRR556174_1.fastq somewhere/SRR556174_2.fastq somewhere/SRR556175_1.fastq somewhere/SRR556175_2.fastq ttm5/read  
ulimit -n 4096
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0
```

FASTQ files should be saved in the read directory in targets. Files compressed with gz, bzip2, and xz formats can be processed.

Results are saved in the ttm2 directory, specified by a.

After the commit of f52457f, tef.pl opens more than 1000 files simultaneously. Setting ulimit -n 4096 may be required before running tef.pl.

## Making a Reference Data Set
### First Run:
When you first run tef.pl, an empty directory specified by the reference name is created.
Save the gzipped FASTA file for the reference into this directory, and then run tef.pl again.
### Second Run:
After the second run of tef.pl, a config file is created in the reference directory.
The config file contains a list of chromosome names in the order they appear in the FASTA file.
If necessary, rename or change the chromosome names in the config file to “NOP”.
For partial contigs or extra chromosomes that are not required, replace the chromosome name with “NOP”.
If there are too many partial contigs to ignore, create a config file with all lines set to “NOP” and then replace the desired chromosome names at the relevant positions.
### Third Run:
On the third run, ped.pl creates the reference data set according to the config file and then proceeds with the main routine.
Once the reference data set is created, it can be reused for additional analyses.
If the reference data is incomplete or corrupted, delete the reference directory and run tef.pl again.

## sliceTE.pl
Usage
```
perl sliceTE.pl path_of_reference_genome head_sequence tail_sequence
```
The sliceTE.pl outputs sequences between head and tail sequences from the FASTA file of the reference genome. Options should be separated by a single space.

The range of TE size can be adjusted on line 54 of sliceTE.pl.

FASTA files compressed by gzip, bzip2, and xz can be processed.

For example,
```
perl sliceTE.pl IRGSP1.0/IRGSP-1.0_genome.fasta.gz TGTTAAATATATATACAAGC AGGTTGCAAGTTAGTTAAGA
```
## sliceRef.pl
The sliceRef.pl outputs sliced sequences from the chromosome sequence in the reference directory.

For example,
```
perl sliceRef.pl IRGSP1.0/chr7 26694795 26698908
```
The chromosome sequence should be the output of tef.pl.

## Update
- 1.6 TEF now outputs read alignments to the genomic sequences corresponding to transposition sites derived from the detected genotypes. 2026-02-10
- 1.5 Now possible to detect transpositions by comparing a single set of read data with a reference genome sequence. 2025-04-23
- 1.4 Improved accuracy of TE detection. 2023-12-14
- 1.3 Insertion direction of TE is shown with same format of version 1.0. 2023-10-09
- 1.2 Correct output transposed positions. Add filter for ambigure TSDs. 2023-09-18
- 1.1 Improve junction method. Large genome data (even human) can be applicable. 2023-08-23
- 1.0 Initial version. 2022-03-29

## Citing TEF
Miyao, A., Yamanouchi, U. Transposable element finder (TEF): finding active transposable elements from next generation sequencing data. BMC Bioinformatics 23, 500 (2022). https://doi.org/10.1186/s12859-022-05011-3

## Lisence
Free of use for academics. For non-academics, license by NARO (same as PED) is required.  
Patent JP 7573862.  

## Author 
Akio Miyao miyao.akio700@naro.go.jp
