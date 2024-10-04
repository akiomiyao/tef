# TEF
## Transposable Element Finder - Detection of active transposable elements from NGS data  
Transposable Element Finder (TEF) is a detection program of active transposable elements (TEs) by comparing two NGS sequences.
TEF does not require TE information to detect transposition.  
Instead, TEF returns transposed transposon ends information and transposed positon on the reference genome.  
Active TEs transposed on the host genome. Because the transposition events are independent, inserted positions of TE should be different between two samples.
Short read sequences from both samples should contain different fusion fragments of genome and TE sequences.
Most TEs make target site duplication (TSD) at inserted position.
TEF detects both ends of TE and inserted position on the genome with information of TSD from NGS short reads.  
TEF runs on Unix (Linux and FreeBSD) machine.  Eight or 16GB memory is enough. One or more TB disk (SSD is better) space is required. To analyze NGS data of human, 4 or more TB disk space is recommended.  

Web page : https://akiomiyao.github.io/tef/  

## Usage

For example,  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0  

Options are separated with comma WITHOUT space.  

The tef.pl is implemented with junction and TSD methods.  
Default is the junction method.  

If tsd_size option is specified, tef.pl will run with the TSD method.  

For example,  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5  

The th option is threshold for filtering same TSD sequence.  
Default = 0.2  
If too many noises are detected, increase th value, *e.g.* th=0.7.  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,th=0.7  

The th option is TSD method specific.  

Because tmp directries in tergets grow huge size, tmp directories will be deleted  
at the end of analysis. If debug=yes is added in options, tmp directories will not  
deleted.  
  
sort_tmp=directory_for_sort is the option for temporary directory for sort command.   
If sort_tmp is specified to fast disk, *e.g.* SSD, sorting will be accelerated.  
If your SSD is not enough size, set the tef directory in HDD space and specify sort_tmp in SDD space.  

For Linux, max process is limited to number of CPU core. For other OS, default process number is 4.  
If you add max_process option, *e.g.* max_process=8, the tef.pl uses 8 cores.  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5,th=0.7,option=clear,sort_tmp=/mnt/ssd/tmp,max_process=8  

Counting function of k-mer is improved. If disk space becomes an issue, add the compress=yes option.  

In addition, if an error is output and the program stops, execute "ulimit -n 4096" before running this program.  
This will increase the limit on the number of files that can be opened at once.  

## Getting tef.pl
```
git clone https://github.com/akiomiyao/tef.git   
```
for updating tef.pl  
```
cd tef  
git pull  
```
  
## Analysis by tef.pl  
At first, directories specified by a and b should be made.  

For example,  
```
cd tef  
mkdir IRGSP1.0  # If tef.pl runs without ref directory, tef.pl makes ref directory and exit.
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

Fastq files should be saved in read directory in targets.  
Files compressed with gz, bzip2, and xz format can be processed.  
  
Results are saved in ttm2 directory, specified by 'a'.  
  
After the commit of f52457f, tef.pl opens more than 1000 files simalteneously.  
Setting  
ulimit -n 4096  
may be required, before running tef.pl  

## Making reference data set
At the first run of tef.pl, the empty directory with specified by reference name is created.  
Save the gz compressed fasta file for the reference into the directory, and then run again tef.pl.  
After the second run of tef.pl, config file is created in the reference directory.  
Config file is the list of chromosome names in order to appearance in fasta file.  
If required, rename and/or change to NOP the chromosome name in the config file.  
If the sequence is partial contig or extra chromosome and is not required, the chromosome name should be replaced to NOP.  
In the case of too many partial contigs which should be ignore, make config file with all lines are NOP and then replace disired chromosome names at related positions.  
At the third run, ped.pl makes reference data set according to the config file and then proceed main routine.  
Once the reference data set is created, data set will be reused for additional analysis.  
If the reference data is imcomplete or broken, delete the reference directory and then run again tef.pl.  

## sliceTE.pl
Usage
perl sliceTE.pl path_of_reference_genome head_sequence tail_sequence

The sliceTE.pl outputs sequences between head and tail sequences from fasta file of reference genome. Options should be separated by single space.  

The range of TE size can be adjusted on line 54 of sliceTE.pl.  

Fasta file compressed by gzip, bzip2 and xz can be processed.  

For example,
```
perl sliceTE.pl IRGSP1.0/IRGSP-1.0_genome.fasta.gz TGTTAAATATATATACAAGC AGGTTGCAAGTTAGTTAAGA
```
## sliceRef.pl
The sliceRef.pl outputs sliced sequence from the chromosome sequence in reference directory.

For example,
```
perl sliceRef.pl IRGSP1.0/chr7 26694795 26698908
```
The chromosome sequence should be the output of tef.pl.

## Update
- 1.4 Improved accuracy of TE detection. 2023-12-14
- 1.3 Insertion direction of TE is shown with same format of version 1.0. 2023-10-09
- 1.2 Correct output transposed positions. Add filter for ambigure TSDs. 2023-09-18
- 1.1 Improve junction method. Large genome data (even human) can be applicable. 2023-08-23
- 1.0 Initial version. 2022-03-29

## Citing TEF
Miyao, A., Yamanouchi, U. Transposable element finder (TEF): finding active transposable elements from next generation sequencing data. BMC Bioinformatics 23, 500 (2022). https://doi.org/10.1186/s12859-022-05011-3

## Lisence
Free of use for academics. For non-academics, license by NARO (same as PED) is required.  
Patent pending JP 2020-217693.  

## Author 
Akio Miyao miyao@affrc.go.jp
