hpg-libs
========

C libraries for the High Performance Genomics (HPG) project:

* aligners: HPC implementation of two widely used sequence aligners: Burrows-Whheler Transform (BWT) and Smith-Waterman (SW)
* bioinfo-libs:  several data format parsers: gff, bam, vcf, bed, ...
* common-libs: several common utilities for logs, string manipulation, workflows execution, ...
* containers: several data structures : array_list, linked_list, ... other have been collected from different sources such as cprops or khash.
* math: Statistics and algebra functions


Each library must be compiled independently using SCons. A static library file (.a) will be generated as a result in c/build/ folder.

BUILDING
Scons is used as building system, to build the library download 'next' branch (for the most recent library), and execute:

 $ scons


GCC4.4+ is the only requirement.


DOCUMENTATION
You can find more info about HPG project and bioinfo-libs at:

 http://wiki.opencb.org/projects/hpg/doku.php


CONTACT
You can contact any of the following developers:
 * Joaquin Tarraga (jtarraga@cipf.es)
 * Cristina Y. Gonzalez (cgonzalez@cipf.es)
 * Ignacio Medina (imedina@cipf.es)
 * Hector Martinez (martineh@cipf.es)
