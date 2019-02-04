# AbundantOTU
Program Name: AbundantOTU+ (AbundantOTU with additional functions)
Version: 0.95b
Released on Feb 1, 2019
Older versions:
AbundantOTU+0.93b was released on Feb 10, 2013
AbundantOTU+0.92b was released on Jan 10, 2012
AbundantOTU+0.91b was released on July 10, 2011
AbundantOTU+0.9b was released on Jan 25, 2011
AbundantOTU2.0 was released on Oct 29, 2010 
AbundantOTU1.0 was released on July 7, 2010
Developer: Yuzhen Ye <yye@indiana.edu>
Affiliation: School of Informatics and Computing, Indiana University, Bloomington

The development of AbundantOTU was supported by NIH grants 1R01HG004908 & 1U01HL098960

AbundantOTU+ is free software under the terms of the GNU General Public License as published by 
the Free Software Foundation.

>> What's new in AbundantOTU+0.95b
   * Add support for paired reads (using -i read1-file -j read2-file)
   * Add scripts (under scripts folder) for processing multiple samples to output an OTU-sample table 
   * Add example files (under msample folder) to demonstrate the steps from multiple samples to OTU-sample table 

>> What's new in AbundantOTU+0.93b
   Bugs fixed--a couple of variables were not initialized properly, which may result in differnet results on different machines of the same dataset.

>> What's new in AbundantOTU+0.92b

   AbundantOTU+0.92b fixed a few bugs

   Added option: -abundantonly
   If you specify -abundantonly, AbundantOTU+ will stop after the consensus alignment step (for inference of abundant species) without attemping the heuristic step to group the rare sequences -- making AbundantOTU+ equivalent to AbundantOTU.
   When to use -abundantonly:
      1) you only care about abundant OTU
      2) you have been using AbundantOTU2.0 and felt it is just right for you
      3) to speedup the calculation--because the last step is skipped 

>> What's new in AbundantOTU+
   AbundantOTU+ works to infer abundant OTUs (using consensus alignment algorithm) and rare OTUs (using heuristic strategy)

>> Notes before you start
  * AbundantOTU+ was tested extensively on linux computers. It should work on linux and unix machines. 
  * How much memory will AbundantOTU+ consume? It uses about 1.5 of the size of the input sequence file (or less). 

>> Introduction

AbundantOTU+ is a tool for fast identification and quantification of abundant Operational Taxonomic Units (OTU). It takes in DNA sequences of 16S rRNA genes (gene fragments) as inputs, and outputs consensus sequences of abundant OTUs. 

AbundantOTU+ (and AbundantOTU) on the web:
  http://omics.informatics.indiana.edu/AbundantOTU+
  (Please check the project home page for updates and newer version of the program)

>> Installation

Simply run install

The executable file "AbundantOTU" will be created under folder bin/

>> Using AbundantOTU+

Usage: type AbundantOTU+ for usages

See an example under examples/, go to this folder, and run,

  ../bin/AbundantOTU+ -i DivergentGSFLXClean.fa -o DivergentGSFLXClean-AbundantOTU+

  Input: DivergentGSFLXClean.fa
  Outputs:
   DivergentGSFLXClean-AbundantOTU.cons      -- a file of consensus sequences, each representing an abundant OTU
   DivergentGSFLXClean-AbundantOTU.clustsize -- a file of the sizes of the OTUs
   DivergentGSFLXClean-AbundantOTU.clust     -- a file of detailed information of the recruitment of the reads

  Note: there is no more *.single files produced by AbundantOTU+ (unless you use -abundantonly option);
        but you may decide the abundant level of the OTUs you are interested by checking the *.clustsize file

>> Using AbundantOTU+ for multiple samples
  Check msample/readme for details; example files can also be found in msample folder.
  Note that you can change parameters when calling otu-sample-table.py to produce the OTU-sample table according
  to your need. 
  Examples, "-m 5" to change minimum read suppport (5) for the OTUs to be included in the table; 
            "-n no" to output read counts (instead of proportions) in the table.

>> Notes
  If AbundantOTU+ claims that there are many sequences with multiple seeds, if means that many of your sequences contain repeats. One possibility is that you forgot to trim your sequences--they contain identical tags at both ends of the sequences.
