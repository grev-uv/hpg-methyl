Direct BWT construction

Kunihiko Sadakane
National Institute of Informatics (NII)
http://researchmap.jp/sada/

Usage:

dbwt filename
  compute the BWT of the file "filename".  Output files are "output.bw" and
  "output.lst" (the rank of the last character).
  The output files can be used for csalib.

Limitation:
  The program will work only if the input file is of size < 4GB.

Acknowledgment:
  The file "sais.c" was written by Yuta Mori.
  http://sites.google.com/site/yuta256/

Algorithm:
  This program is based on the algorithm proposed in
    Daisuke Okanohara, Kunihiko Sadakane. A Linear-Time Burrows-Wheeler Transform 
    Using Induced Sorting. In Proc. of SPIRE, LNCS 5721, pp. 90-101, 2009,
  with some simplification.  For most of inputs, the program uses less than 2.5n bytes
  where n is the length of the input.

http://dx.doi.org/10.1007/978-3-642-03784-9_9

Changes:
2010-07-30: Fixed a bug (the previous one did not work on Ubuntu).
2010-07-17: First release.
