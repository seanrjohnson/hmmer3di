## HMMER3Di - This is a fork of hmmer3.3.2 and easel 0.48 patched to support the Foldseek (3Di) alphabet

![](http://img.shields.io/badge/license-BSD-brightgreen.svg)


__The original hmmer and easel repositories are here__

[github](https://github.com/EddyRivasLab/hmmer). 

and here:
[github](https://github.com/EddyRivasLab/easel).


## Install

```
   % git clone git@github.com:seanrjohnson/hmmer3di.git 
   % autoconf
   % ./configure --prefix /your/install/path
   % make
   % source copy_executables.sh 
   % # copy executables will create a new directory called hmmer3Di, then it will copy
   % # hmmalign, hmmbuild, hmmpress, hmmsearch, hmmscan, and phmmer into that directory
   % # with 3Di_ added to the start of their names. From there you can execute them 
   % # or copy them into your $PATH
``` 

## Background frequencines for 3Di sequences can be found here

[3Di_background_frequencies.txt](3Di_background_frequencies.txt)

## Dirichlet priors calculated from 3Di MSAs
To generate a set of 3Di MSAs, we converted the AlphaFold UniProt Foldseek database (Jumper et al., 2021; van Kempen et al., 2023; Varadi et al., 2022) to a 3Di fasta file. We then looked up every sequence name from the Pfam 35 seed file in the UniProt 3Di fasta file and, for cases where the corresponding sequence was identifiable, extracted the sub-sequence corresponding to the Pfam 35 seed. 3Di seeds from each profile were aligned using MAFFT. MSA columns with more than 10 rows were used to calculate background frequencies and Dirichlet priors using the HMMER3 program esl-mixdchlet fit with options -s 17 9 20. 
[pfam_35_3Di_msa_counts_lb_10.mixdchlet.txt](pfam_35_3Di_msa_counts_lb_10.mixdchlet.txt)

## Changes made to support the 3Di alphabet

A full list of changes can be seen in the following diff:
[https://github.com/seanrjohnson/hmmer3di/compare/2637afc..87a5d15](https://github.com/seanrjohnson/hmmer3di/compare/2637afc..87a5d15)