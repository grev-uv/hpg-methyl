CSA library
http://researchmap.jp/sada/csalib/
http://code.google.com/p/csalib/

Kunihiko Sadakane
National Institute of Informatics (NII)
sada@nii.ac.jp

Usage:

mkcsa filename [-P[n]:{L}:{O}] [-I:{D}:{D2}]
  It constructs a compressed full-text index from the BWT.
  The program reads filename.bw and filename.lst, and outputs
  Psi/BW.  It outputs sampled SA and inverse SA if -I option is given.
  -P[n] specifies the compression method.

    -P1    psi encoded by the gamma function
    -P2    psi encoded by the gamma function and run-length codes
    -P3    BW for DNA
    -P4    BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed.
    -P5    BW using Huffman-shaped wavelet tree.  Bit-vectors are further compressed.
    -P6    psi encoded by the sparse array
    -P7    BW using Huffman-shaped wavelet tree.  Bit-vectors are not compressed (dense array).
    -P10   BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed with sparse array.
    -P11   psi encoded by lengths of zero-runs and one-runs
    -P12   BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed with lengths of zero-runs and one-runs.

  If -I is given, the common data (frequencies of characters, sampled SA) are output.
  Output filename is filename.idx.  This file is compatible with any compression method.

  L, O, D, D2 are options.
     L specifies the sampling rate.
     D specifies the sampling rate of SA
     D2 specifies the sampling rate of inverse SA

  Recomended parameters are
     -P1:256 (good balance of size and speed)
     -P4:512 (good compression)
     -P12:512 (good compression)
     -P3:512 (absolutely good for DNA)

  The input files filename.bw and filename.lst can be created from a file "filename"
  by dbwt or ssss commands.  They are available on the project Web page.

csa filename.idx filename.{ext}
  It tests the index.

unbwt filename.idx filename.{ext}
  It decodes the compressed file.  Output filename is "output.dec".

API:
int csa_read(CSA *SA, int argc, char *argv[]);
  read the index into memory.  The structure SA stores pointers to supported functions.

SA->search(unsigned char *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri)
  search the index for pattern "key" of length "keylen".
  The output is the range of the suffix array [li,ri]
SA->lookup(CSA *csa, rank_t i)
  return SA[i]
SA->inverse(CSA *csa, pos_t suf)
  return ISA[i] (inverse SA)
SA->text(uchar *p,CSA *csa,pos_t i,pos_t j)
  extract T[i,j] to "p"

There are many other functions not listed here.

Compile:
  The library is compiled with gcc4.3 to use the Intel SSE4.2 technology supported
  by Core i7 processor.  SSE4.2 has the popcount instruction, which improves the speed
  of succinct data structures.  The library does not check if the CPU has SSE4.2.
  In this case, the program will terminate with illegal instruction error.
  If the -msse4.2 option is not given, the library will use software popcount.

Changes:
2010-08-10: Fixed a bug in psi1.c.  Improved the speed for successor/predecessor of Psi.
            Added sample programs (left/right maximal frequent patterns).
2010-07-18: Added ruby binding and Pizza&Chile interface.
2010-07-17: First release.
