## HMMER - biological sequence analysis using profile HMMs

[![](https://travis-ci.org/EddyRivasLab/hmmer.svg?branch=develop)](https://travis-ci.org/EddyRivasLab/hmmer)
![](http://img.shields.io/badge/license-BSD-brightgreen.svg)

[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models" (profile HMMs). HMMER is used by many
protein family domain databases and large-scale annotation pipelines,
including [Pfam](http://pfam.xfam.org) and other members of the
[InterPro Consortium](http://www.ebi.ac.uk/interpro/).

__This is a fork of hmmer3.3.2 and easel 0.48 patched to support the foldseek (3di) alphabet__


__The original hmmer and easel repositories are here__
To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).



```
   % git clone git@github.com:seanrjohnson/hmmer3di.git 
   % autoconf
   % ./configure --prefix /your/install/path
   % make
   % make check                 # optional: run automated tests
   % make install               # optional: install HMMER programs, man pages
   % (cd easel; make install)   # optional: install Easel tools
``` 

Executable programs will be installed in `/your/install/path/bin`. If
you leave this optional `./configure` argument off, the default prefix
is `/usr/local`.

Files to read in the source directory:

   * INSTALL - brief installation instructions.
   * Userguide.pdf - the HMMER User's Guide.
 
To get started after installation, see the Tutorial section in the
HMMER User's Guide (Userguide.pdf).


# Guide for adding a new alphabet 

