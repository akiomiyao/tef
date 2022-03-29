# TEF
## Transposable Element Finder - Detection of active transposable elements from NGS data  
Transposable Element Finder (TEF) is a detection program of active transposable elements (TEs) by comparing two NGS sequences.
Active TEs transposed on the host genome. Because the transposition events are independent, inserted positions of TE should be different between two samples.
Short read sequences from both samples should contain different fusion fragments of genome and TE sequences.
Most TEs make target site duplication (TSD) at inserted position.
TEF detects both ends of TE and inserted position on the genome with information of TSD from NGS short reads.

## Usage

For example,  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0  

Options are separated with comma WITHOUT space.  

The tef.pl is implemented with junction method and TSD method.  
Default is the junction method.  

If tsd_size option is specified, tef.pl will run with the TSD method.  

For example,  
perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5  

The th option is threshold for filtering same TSD sequence.  
Default = 0.2  
If too many noises are detected, increase th value, e.g. th=0.7.  

## Getting tef.pl

git clone https://github.com/akiomiyao/tef.git   

for updating tef.pl  
cd tef  
git pull  

## Making reference data set
At the first run of for a reference, the empty directory with specified name is created.  
Save the gz compressed fasta file for the reference into the directory, and then run again tef.pl.  
After the second run of tef.pl, config file is created in the reference directory.  
Config file is the chromosome name list in order to appearance of fasta file.  
If required, rename and/or change to NOP the chromosome name in the config file.  
If the sequence is partial contig or extra chromosome and is not required, the chromosome name should be replaced to NOP.  
At the third run, ped.pl makes reference data set according to the config file and then proceed main routine.  
  
Analysis by tef.pl  
At first, directories specified by a and b should be made.  
  
For example,  
mkdir ttm2  
mkdir ttm2/read  
mkdir ttm5  
mkdir ttm5/read  
  
Fastq files should be saved in read directory in targets.  
Files compressed with gz, bzip2, and xz format can be processed.  
  
Results are saved in ttm2 and ttm5 directory.

## Update
- Initial version 1.0 2022-03-29

## Author 
Akio Miyao miyao@affrc.go.jp
