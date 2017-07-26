/*
	 Copyright 2010, Kunihiko Sadakane, all rights reserved.

	 This software may be used freely for any purpose.
	 No warranty is given regarding the quality of this software.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>
#include "typedef.h"
#include "mman.h"
#include "csa.h"
#include "psi1.h"
#include "lf_dna.h"
#include "lf_wt.h"
#include "psi2.h"
#include "lf_bit.h"


#define display_progressbar(str,i,n) if (i % 10000000 == 0) {fprintf(stderr,"%s %ld/%ld                       \r",str,i/10000000,n/10000000);  fflush(stderr); }

static void csa_error(void)
{
	printf("not supported.\n");
	exit(1);
}

u64 getuint(uchar *s, i64 i, i64 w)
{
	u64 x;
	i64 j;
	s += i*w;
	x = 0;
	for (j=0; j<w; j++) {
		x += ((u64)(*s++)) << (j*8);
	}
	return x;
}

void putuint(uchar *s, i64 i, i64 x, i64 w)
{
	i64 j;
	s += i*w;
	for (j=0; j<w; j++) {
		*s++ = x & 0xff;
		x >>= 8;
	}
}

static int blog(i64 x) // [0,n] の数を格納するには blog(n)+1 ビット必要
{
	int l;
	l = -1;
	while (x>0) {
		x>>=1;
		l++;
	}
	return l;
}

void writeint(int k,i64 x,FILE *f)
{
	int i;
	for (i=k-1; i>=0; i--) {
		fputc(x & 0xff,f);
		x >>= 8;
	}
}

#if 0
i64 readint(int k,FILE *f)
{
	i64 x;
	int i;
	x = 0;
	for (i=0; i<k; i++) {
		x += ((i64)fgetc(f)<<(8*i));
	}
	return x;
}

u64 readuint(int k,FILE *f)
{
	u64 x;
	int i;
	x = 0;
	for (i=0; i<k; i++) {
		x += ((u64)fgetc(f)<<(8*i));
	}
	return x;
}

#endif

void writeuint(int k,u64 x,FILE *f)
{
	int i;
	for (i=k-1; i>=0; i--) {
		fputc(x & 0xff,f);
		x >>= 8;
	}
}

static diskbuf *psi[SIGMA];
void bw_to_psi(FILE *out, CSA *csa, char *fbw, char *flst, int *k)
{
	FILE *in;
	i64 last = 0,i,j;
	i64 C2[SIGMA];
	i64 c;

	in = fopen(flst,"r");
	if (in == NULL) {
		perror("bw_to_psi:");  exit(1);
	}

	fclose(in);

	for (c=0; c<SIGMA; c++) {
		csa->C[c] = 0;
	}

	in = fopen(fbw,"r");
	if (in == NULL) {
		perror("bw_to_psi:");  exit(1);
	}
	csa->n = 0;
	while (1) {
		display_progressbar("reading ",csa->n,0L);
		c = fgetc(in);
		if (c == EOF) break;
		csa->C[c]++;
		csa->n++;
	}
	rewind(in);
	printf("n = %ld last = %ld\n",csa->n,last);

	*k = (blog(csa->n+1)+1+8-1)/8;
	for (c=0; c<SIGMA; c++) {
		psi[c] = open_diskbuf(out,*k);
	}
	for (j=1,c=0; c<SIGMA; c++) {
		C2[c] = j;
		j += csa->C[c];
	}

	for (i = 0; i<=csa->n; i++) {
		display_progressbar("computing psi ",i,csa->n);
		if (i == last) {
			setint_diskbuf(psi[0], 0, i);
		} else {
			c = fgetc(in);
			setint_diskbuf(psi[c], C2[c]++, i);
		}
	}
	fclose(in);
	for (c=0; c<SIGMA; c++) {
		close_diskbuf(psi[c]);
	}
}

i64 write_header(CSA *csa, FILE *f2, int isize)
{
	i64 i,m,k,c;

	m = csa->m;

	k = (blog(csa->n+1)+1+8-1)/8;
	writeint(1,k,f2); /* #bytes of integer */
	writeint(k,csa->n,f2);   /* length of the text */
	//  writeint(k,csa->L,f2); /* interval between two psi values stored explicitly */
	writeint(k,csa->D,f2); /* interval between two SA values stored explicitly */
	writeint(k,csa->D2,f2); /* interval between two inverse SA values stored explicitly */
	writeint(k,SIGMA,f2);   /* alphabet size */
	writeint(k,csa->m,f2);   /* the number of distinct characters in the text  */
	isize += 1+5*k;

	for (i = 0; i < m; i++) {
		c = csa->AtoC[i];
		writeint(1,c,f2); /* characters appeared in the text */
		writeint(k,csa->C[c],f2); /* frequency of characters */
	}
	isize += m*(k+1);

	return isize;
}

i64 write_sa(CSA *csa, FILE *f2, int isize)
{
	i64 i,n,D,k;
	n = csa->n;  D = csa->D;  k = csa->k;
	printf("write_sa n=%ld D=%ld k=%ld\n",n,D,k);
	writeint(1,k,f2); /* #bytes of integer */
	writeint(k,D,f2); /* interval between two SA values stored explicitly */
	isize += 1+k;
	for (i = 0; i <= n/D; i++) {
		display_progressbar("writing sa ",i,n/D);
		//    printf("sa[%ld] = %ld\n",i,getuint(sa->SA,i,k));
		writeint(k,getuint(csa->SA,i,k),f2);
		isize += k;
	}
	return isize;
}

i64 write_isa(CSA *csa, FILE *f2, int isize)
{
	i64 i,n,D2,k;
	n = csa->n;  D2 = csa->D2;  k = csa->k;
	writeint(1,k,f2); /* #bytes of integer */
	writeint(k,D2,f2); /* interval between two SA values stored explicitly */
	isize += 1+k;
	for (i = 0; i <= (n-1)/D2; i++) {
		display_progressbar("writing isa ",i,(n-1)/D2);
		writeint(k,getuint(csa->ISA,i,k),f2);
		isize += k;
	}
	return isize;
}

static void read_header(CSA *csa, uchar **map)
{
	i64 i,k,m,c;
	uchar *p;

	//  printf("read_header\n");
	p = *map;

	k = csa->k = getuint(p,0,1);  p += 1;
	//  printf("k = %ld\n",k);
	csa->n = getuint(p,0,k);  p += k; /* length of the text */
	csa->D = getuint(p,0,k);  p += k; /* interval between two SA values stored explicitly */
	csa->D2 = getuint(p,0,k);  p += k; /* interval between two inverse SA values stored explicitly */
	m = getuint(p,0,k);  p += k;
	if (m != SIGMA) {   /* alphabet size */
		printf("error sigma=%ld\n",m);
	}
	csa->m = m = getuint(p,0,k);  p += k; /* the number of distinct characters in the text */

	for (i = 0; i < SIGMA; i++) csa->C[i] = 0;
	for (i = 0; i < m; i++) {
		c = getuint(p,0,1);  p += 1;
		csa->C[c] = getuint(p,0,k);  p += k;
	}

	*map = p;
}

static void read_sa(CSA *csa, uchar **map)
{
	i64 k,n;
	uchar *p;

	p = *map;

	n = csa->n;
	k = getuint(p,0,1);  p += 1;
	csa->D = getuint(p,0,k);  p += k;

	printf("read_sa n=%ld k=%ld D=%d\n",n,k,csa->D);

	csa->SA = p;
	p += (n / csa->D+1)*k;

	*map = p;
}

static void read_isa(CSA *csa, uchar **map)
{
	i64 k,n;
	uchar *p;

	p = *map;

	n = csa->n;
	k = getuint(p,0,1);  p += 1;
	csa->D2 = getuint(p,0,k);  p += k;

	printf("read_isa n=%ld k=%ld D2=%d\n",n,k,csa->D2);

	csa->ISA = p;
	p += ((n-1) / csa->D2+1)*k;

	*map = p;
}

void csa_options(CSA *csa, char *p)
{
	int idx_id;
	csa->D = 128;
	csa->D2 = csa->D*2;
	if (p[0] == 0) return;
	if (p[0] != ':') {
		sscanf(p,"%d",&idx_id);
		printf("idx_id = %d\n",idx_id);
		while (p[0] != ':' && p[0] != 0) p++;
	}
	if (p[0] == 0) return;
	p++;
	if (p[0] != ':') {
		sscanf(p,"%d",&csa->D);
		printf("D = %d\n",csa->D);
		csa->D2 = csa->D*2;
		while (p[0] != ':' && p[0] != 0) p++;
	}
	if (p[0] == 0) return;
	p++;
	if (p[0] != ':') {
		sscanf(p,"%d",&csa->D2);
		printf("D2 = %d\n",csa->D2);
		while (p[0] != ':' && p[0] != 0) p++;
	}
	if (p[0] == 0) return;
	p++;
}

void psi_options(CSA *csa, char *p)
{
	csa->id = ID_DIFF_GAMMA;
	//  if (p[0] == 0) return;
	if (p[0] != 0 && p[0] != ':') {
		sscanf(p,"%d",&csa->id);
		while (p[0] != ':' && p[0] != 0) p++;
	}
	switch (csa->id) {
		case ID_DIFF_GAMMA:
		case ID_DIFF_GAMMA_RL:
		case ID_DIFF_GAMMA_SPARSE:
		case ID_DIFF_GAMMA_RL_SPARSE:
		case ID_DIFF_GAMMA_RR:
			psi1_options(csa, p);
			break;
		case ID_BWT_DNA:
			lf_dna_options(csa, p);
			break;
		case ID_BWT_BIT:
			lf_bit_options(csa, p);
			break;
		case ID_BWT_WT:
		case ID_BWT_WT_HUF:
		case ID_BWT_WT_DENSE:
		case ID_BWT_WT_SPARSE4:
		case ID_BWT_WT_RR:
			lf_wt_options(csa, p);
			break;
#if 0
		case ID_BWT_HUF:
			lf_bwt_options(csa, p);
			break;
#endif
		case ID_SPARSE4:
			psi2_options(csa, p);
			break;
		default:
			printf("error unknown psi_id %d\n",csa->id);
			exit(1);
			break;
	}
}

void csa_new_from_bwt(CSA csa, char *fname, char *fidx, int psi_id, int idx_id, bool coded)
{
	int k;
	i64 i,j,v,m;
	FILE *f2;
	i64 psize,isize;
	i64 n;

	psi_id = csa.id;
	if (psi_id >= 0) {
		printf("create psi: id=%d\n",psi_id);
	}
	if (idx_id >= 0) {
		printf("create idx: id=%d D=%d D2=%d\n",idx_id,csa.D,csa.D2);
	}

	psize = 0;

	if (psi_id >= 0) {
		switch (psi_id & 0x3f) {
			case ID_DIFF_GAMMA:
			case ID_DIFF_GAMMA_RL:
			case ID_DIFF_GAMMA_SPARSE:
			case ID_DIFF_GAMMA_RL_SPARSE:
				psize = psi1_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("Psi   %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
			case ID_DIFF_GAMMA_RR:
				psize = psi12_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("Psi   %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
			case ID_BWT_DNA:
				psize = lf_dna_makeindex(&csa, fname, coded);
				printf("n     %ld\n",csa.n);
				printf("BW    %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
			case ID_BWT_BIT:
				psize = lf_bit_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("BW    %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
			case ID_BWT_WT:
			case ID_BWT_WT_HUF:
			case ID_BWT_WT_DENSE:
			case ID_BWT_WT_SPARSE4:
			case ID_BWT_WT_RR:
				psize = lf_wt_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("BW    %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
#if 0
			case ID_BWT_HUF:
				psize = lf_bwt_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("BW    %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
#endif
			case ID_SPARSE4:
				psize = psi2_makeindex(&csa, fname);
				printf("n     %ld\n",csa.n);
				printf("Psi   %ld bytes (%1.3f bpc)\n",
						psize,(double)psize*8/csa.n);
				break;
			default:
				printf("psi_id = %d\n",psi_id);
				exit(1);
		}
	}

	csa.k = (blog(csa.n+1)+1+8-1)/8;

	for (i=0; i<SIGMA; i++) csa.CtoA[i] = -1;
	//  csa.K[-1+1] = 0;
	csa.K[-1+1] = 1;
	for (m=0,v=1,i=0; i<SIGMA; i++) {
		if (csa.C[i]>0) {
			csa.AtoC[m] = i;
			csa.CtoA[i] = m;
			csa.K[m+1] = v;
			//      printf("i=%ld v = %ld C[i] = %ld\n",i,v,csa.C[i]);
			v += csa.C[i];
			m++;
		}
	}
	csa.K[m+1] = v;
	csa.m = m;

	if (csa.D >= csa.n) {
		printf("D=%d >= n=%ld\n",csa.D,csa.n);
		exit(0);
	}
	if (csa.D2 >= csa.n) {
		printf("D2=%d >= n=%ld\n",csa.D2,csa.n);
		exit(0);
	}

	if (idx_id >= 0) {
		n = csa.n;
		k = csa.k;
		////  compute SA and ISA
		if (csa.D > 0) csa.SA = (uchar *) mymalloc(((n-1)/csa.D+1+1)*k);
		if (csa.D2 > 0) csa.ISA = (uchar *) mymalloc(((n-1)/csa.D2+1+1)*k);
		if (csa.D == 0 && csa.D2 == 0) goto brk;

		switch (psi_id & 0x3f) {
			case ID_DIFF_GAMMA:
			case ID_DIFF_GAMMA_RL:
			case ID_DIFF_GAMMA_SPARSE:
			case ID_DIFF_GAMMA_RL_SPARSE:
			case ID_SPARSE4:
			case ID_DIFF_GAMMA_RR:
				j = 0;
				for (i=0; i<=n; i++) {
					display_progressbar("making sa ",i,n);
					j = csa.psi(&csa,j);
					//  sa[j] = i;
					if (csa.D > 0 && j % csa.D == 0) {
						putuint(csa.SA,j / csa.D,i,k);
					}
					if (csa.D2 > 0 && i % csa.D2 == 0) {
						putuint(csa.ISA,i / csa.D2,j,k);
					}
				}
				//      putuint(csa.SA,0,n,k);
				break;
			case ID_BWT_DNA:
			case ID_BWT_BIT:
			case ID_BWT_WT:
			case ID_BWT_WT_HUF:
			case ID_BWT_WT_DENSE:
			case ID_BWT_WT_SPARSE4:
			case ID_BWT_WT_RR:
			case ID_BWT_HUF:
				j = 0;
				for (i=n-1; i>=0; i--) {
					display_progressbar("making sa ",i,n);
					v = csa.LF(&csa,j);
					//        printf("LF[%ld] = %ld\n",j,v);
					j = v;
					if (csa.D > 0 && j % csa.D == 0) putuint(csa.SA, j/csa.D , i, k);
					if (csa.D2 > 0 && i % csa.D2 == 0) putuint(csa.ISA, i/csa.D2, j, k);
				}
				putuint(csa.SA,0,n,k);
				break;
			default:
				break;
		}
brk:
		////      write idx
		f2 = fopen(fidx,"wb"); /* directory */
		if (f2 == NULL) {
			perror("csa2_new1: ");
			exit(1);
		}

		isize = 0;

		writeint(4,VERSION,f2); /* version */
		isize += 4;

		writeint(1,ID_HEADER,f2); // header ID
		isize += 1;
		isize = write_header(&csa, f2, isize);

		if (csa.D > 0) {
			writeint(1,ID_SA,f2);
			isize += 1;
			isize = write_sa(&csa, f2, isize);
		}

		if (csa.D2 > 0) {
			writeint(1,ID_ISA,f2);
			isize += 1;
			isize = write_isa(&csa, f2, isize);
		}


		fclose(f2);

		if (csa.D > 0) free(csa.SA);
		if (csa.D2 > 0) free(csa.ISA);

		printf("Total %ld bytes (%1.3f bpc)\n",(psize+isize),
				(double)(psize+isize)*8/csa.n);
	}
	free(fidx);
}

void csa_new_from_bwt_gnu_bwt_wrapper(const char *directory, const char *name) {
	i64 i;
	int psi_id, idx_id;
	char *fname, *fidx;
	CSA csa;

	for (i=0; i<SIGMA+2; i++) csa.C[i] = 0;

	idx_id = 0;
	psi_id = 0;

	char csa_opt[10];
	strcpy(csa_opt, ":64:64");
	char psi_opt[10];
	strcpy(psi_opt, "3:64");

	csa_options(&csa, csa_opt);
	psi_options(&csa, psi_opt);

	fname = (char *) malloc(500 * sizeof(char));
	fidx = (char *) malloc(500 * sizeof(char));
	fname[0]='\0';
	fidx[0]='\0';

	if (name==NULL) {
		strcat(fname, directory);
		strcat(fname, "/output");
  } else {
		strcat(fname, directory);
		strcat(fname, "/");
		strcat(fname, name);
	}

	strcat(fidx, fname);
	strcat(fidx, ".idx");

	csa_new_from_bwt(csa, fname, fidx, psi_id, idx_id, true);
}

void csa_new_from_bwt_wrapper(int argc, char *argv[]) {
	i64 i;
	char *p;
	int k;
	int psi_id, idx_id;
	char *fname,*fidx = 0;
	CSA csa;

	for (i=0; i<SIGMA+2; i++) csa.C[i] = 0;

	fname = NULL;
	psi_id = idx_id = -1;

	for (i=1; i<argc; i++) {
		p = argv[i];
		if (p[0] == '-') {
			p++;
			switch (toupper(p[0])) {
				case 'I':
					// -I[n]:[D]:[D2]
					p++;
					idx_id = 0;
					csa_options(&csa, p);
					break;
				case 'P':
					// -P[n]:[L]
					p++;
					psi_id = 0;
					psi_options(&csa, p);
					break;
				default:
					printf("??? %s\n",argv[i]);
					exit(1);
			}
		} else {
			fname = argv[i];
			k = strlen(fname);
			fidx = (char *) mymalloc(k+5);
			sprintf(fidx,"%s.idx",fname);
		}
	}

	if (fname == NULL) {
		printf("no input file.\n");
		exit(0);
	}

	csa_new_from_bwt(csa, fname, fidx, psi_id, idx_id, false);

}

i64 read_idx(CSA *csa, char *fidx)
{
	i64 isize;
	i64 i,id;
	uchar *p,*q;
	MMAP *mapidx;

	//  printf("read_idx: %s\n",fidx);

	mapidx = mymmap(fidx);
	csa->mapidx = (void *)mapidx;
	p = q = (uchar *) mapidx->addr;
	if (p == NULL) {
		perror("mmap2\n");
		exit(1);
	}
	isize = mapidx->len;

	i = getuint(p,0,4);  p += 4; /* version */
	if (i != VERSION) {
		printf("read_csa: Version %ld is not supported.\n",i);
		exit(1);
	}

	while (1) {
		if (p - q >= isize) break;
		id = getuint(p,0,1);  p += 1;
		//printf("header ID=%ld\n",id);
		switch (id) {
			case ID_HEADER:
				read_header(csa, &p);
				break;

			case ID_SA:
				read_sa(csa, &p);
				break;

			case ID_ISA:
				read_isa(csa, &p);
				break;

			default:
				printf("read_csa: ID %ld is not supported.\n",id);
				break;
		}
	}

	return isize;
}

int csa_read(CSA *csa, int argc, char *argv[])
{
	i64 i,m;
	i64 psize,isize;
	i64 k;
	psi1 *ps;
	psi2 *ps2;
	lf_dna *lf_d;
	lf_wt *lf_w;
	char *fname;
	int ac;

	csa->psi = (rank_t (*)(CSA *, rank_t))csa_error;
	csa->LF = (rank_t (*)(CSA *,rank_t))csa_error;
	csa->rankc = (rank_t (*)(CSA *,rank_t, int))csa_error;
	csa->BW = (int (*)(CSA *,rank_t))csa_error;

	csa->BW_rank = csa_BW_rank;
	csa->lookup = csa_lookup;
	csa->inverse = csa_inverse;
	csa->text = csa_text;
	csa->substring = csa_substring;
	csa->substring_lf = csa_substring_lf_naive;
	csa->T = csa_T;
	csa->head = csa_head_rank;
	csa->search = csa_search;
	csa->searchsub = csa_searchsub;
	csa->child_l = csa_child_l;
	csa->child_r = csa_child_r;

	isize = 0;  psize = 0;

	for (ac = 0; ac < argc; ac++) {

		fname = argv[ac];
		//    printf("argv[%d] = %s\n",ac,argv[ac]);
		k = strlen(fname);

		printf("%s\n", fname+k-4);

		if (strcmp(fname+k-4,".idx") == 0) {	
			isize = read_idx(csa,fname);
		}
		if (strcmp(fname+k-4,".bwd") == 0) {
			////       read bw
			lf_dna_read(csa, fname);
			lf_d = (lf_dna *) csa->psi_struc;
			psize = lf_d->psize;
			//      printf("psize %ld\n",psize);
		}
		if (strcmp(fname+k-4,".bwb") == 0) {
			////       read bw
			lf_bit_read(csa, fname);
			lf_d = (lf_dna *) csa->psi_struc;
			psize = lf_d->psize;
			//      printf("psize %ld\n",psize);
		}
		////
		if (strcmp(fname+k-4,".wtd") == 0
				|| strcmp(fname+k-4,".whd") == 0
				|| strcmp(fname+k-4,".wda") == 0
				|| strcmp(fname+k-4,".wxd") == 0
				|| strcmp(fname+k-4,".wsa") == 0) {
			////       read bw
			lf_wt_read(csa, fname);
			lf_w = (lf_wt *) csa->psi_struc;
			psize = lf_w->psize;
			//      printf("psize %ld\n",psize);
		}
		////
#if 0
		if (strcmp(fname+k-4,".bhd") == 0) {
			////       read bw
			lf_bwt_read(csa, fname);
			lf_w = csa->psi_struc;
			psize = lf_w->psize;
			//      printf("psize %ld\n",psize);
		}
#endif
		////
		if (strcmp(fname+k-4,".psd") == 0
				|| strcmp(fname+k-4,".prd") == 0
				|| strcmp(fname+k-4,".prs") == 0
				|| strcmp(fname+k-4,".pss") == 0) {
			////       read psi  
			psi1_read(csa, fname);
			ps = (psi1 *) csa->psi_struc;
			psize = ps->psize;
			//      printf("psize %ld\n",psize);
			////
		}
		if (strcmp(fname+k-4,".pxd") == 0) {
			////       read psi  
			psi1_read(csa, fname);
			ps = (psi1 *) csa->psi_struc;
			psize = ps->psize;
			//      printf("psize %ld\n",psize);
			////
		}
		if (strcmp(fname+k-4,".psa") == 0) {
			////       read psi  
			psi2_read(csa, fname);
			ps2 = (psi2 *) csa->psi_struc;
			psize = ps2->psize;
			//      printf("psize %ld\n",psize);
			////
		}

	}

	if (csa->D == 0) csa->lookup = (rank_t (*)(CSA *, rank_t))csa_error;
	if (csa->D2 == 0) csa->inverse = (rank_t (*)(CSA *, rank_t))csa_error;

	for (i=0; i<SIGMA; i++) csa->CtoA[i] = -1;
	//  csa->K[-1+1] = 0;
	csa->K[-1+1] = 1;
	for (m=0,k=1,i=0; i<SIGMA; i++) {
		if (csa->C[i]>0) {
			csa->AtoC[m] = i;
			csa->CtoA[i] = m;
			csa->K[m+1] = k;
			//      printf("i=%ld K[%ld] = %ld C[i] = %ld\n",i,m+1,k,csa->C[i]);
			k += csa->C[i];
			m++;
		}
	}
	csa->K[m+1] = k;

	csa->psize = psize;  csa->isize = isize;

#if 0
	fprintf(stderr,"D=%d (stores SA for every D)\n",csa->D);
	//  fprintf(stderr,"L=%d (directory for Psi)\n",csa->l);
	//  fprintf(stderr,"L=%d (directory for Psi)\n",csa->L);
	fprintf(stderr,"n     %ld\n",csa->n);
	fprintf(stderr,"Psi   %ld bytes (%1.3f bpc)\n",
			psize,(double)psize*8/csa->n);
	fprintf(stderr,"Total %ld bytes (%1.3f bpc)\n",psize+isize,
			(double)(psize+isize)*8/csa->n);
#endif

	return 0;
}

/* computes the number of characters smaller than T[SA[i]] */
int csa_head_rank(CSA *csa,rank_t i)
{
	i64 l,r,m;
#ifdef DEBUG
	if (i > csa->n || i < 1) {
		printf("error psi_get i=%u n=%u\n",i,csa->n);
		exit(1);
	}
#endif
	if (i == 0) return -1; // EOF
	l = 0; r = csa->m-1;
	while (l <= r) {
		m = (l + r) / 2;
		if (csa->K[m+1] <= i) l = m + 1; else r = m - 1;
	}
	return l - 1;
}

///////////////////////////////////////////
//  int csa_T(CSA *csa,rank_t i)
//    returns T[SA[i]] (head character of suffix i)
///////////////////////////////////////////
int csa_T(CSA *csa,rank_t i)
{
	int c;
	if (i < 1 || i > csa->n) return -1;
	/* computes the number of characters smaller than T[SA[i]] */
	c = csa->head(csa,i);
	//  if (c == -1) return -1;
	/* converts the number to character code */
	return csa->AtoC[c];
}

///////////////////////////////////////////
//  i64 csa_substring(uchar *p,CSA *csa,rank_t r,i64 len)
//    extracts substring [s,s+len-1] into the buffer p using psi
//    and returns its length
//    length is shorter than len if (s+len-1 > n)
//    where s = SA[r]
///////////////////////////////////////////
i64 csa_substring(uchar *p,CSA *csa,rank_t r,i64 len)
{
	i64 i;
	i = 0;
	while (i < len) {
		if (r == 0) {
			*p++ = 0; // $
			break;
		}
		*p++ = csa->T(csa,r);
		r = csa->psi(csa,r);
		i++;
	}
	return i;
}

///////////////////////////////////////////
//  i64 csa_substring_lf(uchar *p,CSA *csa,rank_t r,i64 len)
//    extracts substring [s-len,s-1] into the buffer p[-len,-1] using BW
//    where s = SA[r]
//  NOTE: the parameters are not compatible with csa_substring!!!
///////////////////////////////////////////
i64 csa_substring_lf(uchar *p,CSA *csa,rank_t r,i64 len)
{
	i64 i,r2;
	int c;

	p += len; // !!!

	i = 0;
	while (i < len) {
		c = csa->BW_rank(csa,r,&r2);
		if (c == -1) {
			*--p = 0; // $
			break;
		}
		*--p = csa->AtoC[c];
		r = csa->K[c+1]-1 + r2;
		i++;
	}
	return i;
}

i64 csa_substring_lf_naive(uchar *p,CSA *csa,rank_t r,i64 len)
{
	i64 i;
	int c;
	i = 0;
	while (i < len) {
		r = csa_BW_LF_by_psi(csa,r, &c);
		if (c == -1) {
			*--p = 0; // $
			break;
		}
		*--p = csa->AtoC[c];
		i++;
	}
	return i;
}

///////////////////////////////////////////
//  rank_t csa_psi_pred_naive(CSA *csa, rank_t pr, rank_t l, rank_t r)
//    returns the index of the predessor of pr in Psi[l,r] using binary search
//            (largest x in [l,r] such that Psi[x] <= pr,
//             l-1 if such x does not exist.)
///////////////////////////////////////////
rank_t csa_psi_pred_naive(CSA *csa, rank_t pr, rank_t l, rank_t r)
{
	i64 m,tmp,ans;
	//i64 l0 = l, r0 = r;
	//  if (pr < csa->psi(csa,l) || pr > csa->psi(csa,r)) {
	//    printf("pred_naive: pr = %ld psi[%ld] = %ld psi[%ld] = %ld\n",
	//           pr,l,csa->psi(csa,l),r,csa->psi(csa,r));
	//    exit(1);
	//  }
	ans = l-1;
	//  printf("pred_naive(%ld, [%ld, %ld])\n",pr, l, r);
	while (l <= r) {
		m = (l+r) / 2;
		tmp = csa->psi(csa,m);
		//    printf("l = %ld r = %ld m = %ld psi = %ld\n", l,r,m, tmp);
		if (tmp == pr) {
			ans = m;
			break;
		}
		if (tmp < pr) {
			ans = m;
			l = m+1;
		} else {
			//      ans = m;
			r = m-1;
		}
	}
//end:
	//  if (pr < csa->psi(csa,l0) || pr > csa->psi(csa,r0)) {
	//  }
	//  if (ans < l0) {
	//    printf("????????????????????????\n");
	//    ans = l0;
	//  }
	//  if (ans > r0) {
	//    printf("????????????????????????\n");
	//    ans = r0;
	//  }
	//  printf("pred(%ld) = %ld\n", pr, ans);
	return ans;
}

///////////////////////////////////////////
//  rank_t csa_psi_succ_naive(CSA *csa, rank_t pl, rank_t l, rank_t r)
//    returns the index of the successor of pl in Psi[l,r] using binary search
//            (smallest x in [l,r] such that Psi[x] >= pl,
//             l-1 if such x does not exist.)
///////////////////////////////////////////
rank_t csa_psi_succ_naive(CSA *csa, rank_t pl, rank_t l, rank_t r)
{
	i64 m,tmp,ans;
	//i64 l0 = l, r0 = r;
	//  if (pl < csa->psi(csa,l) || pl > csa->psi(csa,r)) {
	//    printf("succ_naive: pl = %ld psi[%ld] = %ld psi[%ld] = %ld\n",
	//           pl,l,csa->psi(csa,l),r,csa->psi(csa,r));
	//    exit(1);
	//  }
	ans = r+1;
	//  printf("succ_naive(%ld, [%ld, %ld])\n",pl, l, r);
	while (l <= r) {
		m = (l+r) / 2;
		tmp = csa->psi(csa,m);
		//printf("l = %ld r = %ld m = %ld psi = %ld\n", l,r,m, tmp);
		if (tmp == pl) {
			ans = m;
			break;
		}
		if (tmp > pl) {
			ans = m;
			r = m-1;
		} else {
			l = m+1;
		}
	}
	//  if (pl < csa->psi(csa,l0) || pl > csa->psi(csa,r0)) {
	//    printf("succ(%ld) = %ld\n", pl, r+1);
	//  }
	//  if (ans < l0) {
	//    printf("????????????????????????\n");
	//    ans = l0;
	//  }
	//  if (ans > r0) {
	//    printf("????????????????????????\n");
	//    ans = r0;
	//  }
	//  printf("succ(%ld) = %ld\n", pl, ans);
	return ans;
}



///////////////////////////////////////////
//  rank_t csa_inverse(CSA *csa, pos_t suf)
//    returns SA^{-1}[i] using Psi
///////////////////////////////////////////
rank_t csa_inverse(CSA *csa, pos_t suf)
{
	pos_t p;
	rank_t r;
	i64 D2;

	D2 = csa->D2;

	p = (suf/D2)*D2+1;
	r = getuint(csa->ISA,suf/D2,csa->k);

	while (p < suf+1) {
		r = csa->psi(csa,r);
		p++;
	}
	return r;
}

///////////////////////////////////////////
//  rank_t csa_inverse_lf(CSA *csa, pos_t suf)
//    returns SA^{-1}[i] using BW
///////////////////////////////////////////
rank_t csa_inverse_lf(CSA *csa, pos_t suf)
{
	i64 p,pos;
	i64 t;

	t = csa->D2 - 1;

	if (suf+t < csa->n) {
		//    p = (suf+t) & (~t);
		p = ((suf+t) / (t+1)) * (t+1);
		pos = getuint(csa->ISA,p / (t+1),csa->k);
		while (p > suf) {
			pos = csa->LF(csa,pos);
			p--;
		}
	} else {
		p = csa->n;
		pos = 0;
		while (p > suf) {
			pos = csa->LF(csa,pos);
			p--;
		}
	}
	return pos;
}

/* calculate SA[i] */
///////////////////////////////////////////
//  pos_t csa_lookupf(CSA *csa, rank_t i)
//    returns SA[i] using Psi
///////////////////////////////////////////
pos_t csa_lookup(CSA *csa, rank_t i)
{
	i64 v,D;
	v = 0;  D = csa->D;
	while (i % D !=0) {
		i = csa->psi(csa,i);
		v++;
	}
	i = i / D;
	return getuint(csa->SA,i,csa->k)-v;
}

/* calculate SA[i] */
///////////////////////////////////////////
//  pos_t csa_lookup_lf(CSA *csa, rank_t i)
//    returns SA[i] using BW
///////////////////////////////////////////
pos_t csa_lookup_lf(CSA *csa, rank_t i)
{
	i64 v, D;
	v = 0;  D = csa->D;
	while (i % D !=0) {
		i = csa->LF(csa,i);
		v++;
	}
	i = i / D;
	return (getuint(csa->SA,i,csa->k) + v) % (csa->n+1);
}

///////////////////////////////////////////
//  void csa_text(uchar *p,CSA *csa, pos_t s, pos_t t)
//    extracts substring [s,t] into the buffer p
///////////////////////////////////////////
void csa_text(uchar *p, CSA *csa, pos_t i, pos_t j)
{
	pos_t pos;
	i64 len = j-i+1;
	pos = csa->inverse(csa,i);
	csa->substring(p,csa,pos,len);
}

///////////////////////////////////////////
//  void csa_text_lf(uchar *p,CSA *csa, pos_t s, pos_t t)
//    extracts substring [s,t] into the buffer p using BW
///////////////////////////////////////////
void csa_text_lf(uchar *p,CSA *csa, pos_t s, pos_t t)
{
	i64 i;
	i64 j,c;

#if DEBUG
	if (s < 1 || t > csa->n || s > t) {
		printf("csa_text: out of range [%d,%d] n=%d\n",s,t,csa->n);
	}
#endif
	//  lf = (lf_dna *)csa->psi_struc;

	i = csa->inverse(csa,(t+1));
	for (j=t-s; j>=0; j--) {
		//    printf("rank %ld\n",i);
		c = csa->BW(csa,i);
		if (c == -1) {
			printf("csa_text_lf: ??? j = %ld i = %ld c = %ld  s=%ld t=%ld \n",
					j,i,c,s,t);
			exit(1);
		}
		p[j] = csa->AtoC[c];
		i = csa->LF(csa,i);
	}
}

///////////////////////////////////////////
//  rank_t csa_selectc(CSA *csa, rank_t i, int c)
//    returns select_c(i) using rankc
///////////////////////////////////////////
rank_t csa_selectc(CSA *csa, rank_t i, int c)
{
	i64 l,r,m;

	l = 0; r = csa->n;
	while (l <= r) {
		m = (l+r)/2;
		if (i <= csa->rankc(csa,m,c)) r = m-1;  else l = m+1;
	}
	return l;
}

///////////////////////////////////////////
//  rank_t csa_psi_by_rankc_naive(CSA *csa, rank_t i)
//    returns Psi[i] using BW
///////////////////////////////////////////
rank_t csa_psi_by_rankc_naive(CSA *csa, rank_t i)  /* 0 <= i <= n */
{
	int c;

	if (i==0) return getuint(csa->ISA,0,csa->k); // last

	c = csa_head_rank(csa,i);
	i -= csa->K[c+1];
	i++;

	return csa->selectc(csa,i,c);
}

///////////////////////////////////////////
//  int csa_BW_rank(CSA *csa,i64 i, rank_t *r)
//    returns c = BW[i]
//    *r = rank_c(i)
///////////////////////////////////////////
int csa_BW_rank(CSA *csa,i64 i, rank_t *r)
{
	int c;
	c = csa->BW(csa, i);
	if (c == -1) {
		if (r != NULL) *r = -1;
	} else {
		if (r != NULL) *r = csa->rankc(csa,i,c);
	}
	return c;
}

///////////////////////////////////////////
//  rank_t csa_LF_by_psi(CSA *csa, rank_t i)
//    returns LF[i] using Psi
//  This function is slow (O(sigma log n) time)
///////////////////////////////////////////
rank_t csa_LF_by_psi(CSA *csa, rank_t i)
{
	int c;
	rank_t l,r,x = 0;
	for (c = 0; c < csa->m; c++) {
		l = csa->K[c+1];  r = csa->K[c+2]-1;
		x = csa->psi_pred(csa,i,l,r);
		if (csa->psi(csa,x) == i) break;
	}
	if (c == csa->m) {
		c = -1;
		x = 0; // LF[last] = 0;
	}
	return x;
}

///////////////////////////////////////////
//  rank_t csa_LF(CSA *csa, rank_t i)
//    returns LF[i] using BW
///////////////////////////////////////////
rank_t csa_LF(CSA *csa, rank_t i)  /* 0 <= i <= n */
{
	int c;
	i64 j;

	c = csa->BW(csa,i);
	if (c == -1) return 0; // last
	//  printf("csa_LF: i=%ld c=%d K=%ld\n",i,c,csa->K[c+1]);
	j = csa->K[c+1]-1 + csa->rankc(csa, i, c);

	return j;
}

///////////////////////////////////////////
//  rank_t csa_BW_LF_by_psi(CSA *csa, rank_t i, int *cc)
//    returns LF[i] using Psi
//    *cc = BW[i]
//  This function is slow (O(sigma log n) time)
///////////////////////////////////////////
rank_t csa_BW_LF_by_psi(CSA *csa, rank_t i, int *cc)
{
	int c;
	rank_t l,r,x = 0;
	for (c = 0; c < csa->m; c++) {
		l = csa->K[c+1];  r = csa->K[c+2]-1;
		x = csa->psi_pred(csa,i,l,r);
		if (csa->psi(csa,x) == i) break;
	}
	if (c == csa->m) {
		c = -1;
		x = 0; // LF[last] = 0;
	}
	*cc = c;
	return x;
}

///////////////////////////////////////////
//  i64 csa_searchsub_lf(int c, CSA *csa, rank_t *ll, rank_t *rr)
//    returns 0 if c:key exist, -1 otherwise
//      where key is the substring of length keylen represented by range [*li, *ri]
//    returned [*li, *ri] is the range of c:key
///////////////////////////////////////////
rank_t csa_searchsub_lf(int c, CSA *csa, rank_t *ll, rank_t *rr)
{
	i64 l,r;

	c = csa->CtoA[c];
	if (c == -1) {
		*ll = 1;  *rr = 0;
		return -1;
	}
	l = *ll;  r = *rr;
	l = csa->K[c+1]-1 + csa->rankc(csa, l-1,c) + 1;
	r = csa->K[c+1]-1 + csa->rankc(csa, r,c);
	if (l > r) {
		*ll = 1;  *rr = 0;
		return -1;
	}

	*ll = l;  *rr = r;
	return 0;
}

///////////////////////////////////////////
//  i64 csa_searchsub(int c, CSA *csa, rank_t *ll, rank_t *rr)
//    returns 0 if c:key exist, -1 otherwise
//      where key is the substring of length keylen represented by range [*li, *ri]
//    returned [*li, *ri] is the range of c:key
///////////////////////////////////////////
i64 csa_searchsub(int c, CSA *csa, rank_t *ll, rank_t *rr)
{
	i64 l,r,l0,r0;

	c = csa->CtoA[c];
	if (c == -1) {
		*ll = 1;  *rr = 0;
		return -1;
	}
	l = *ll;  r = *rr;
	r0 = csa->K[c+2]-1;  l0 = csa->K[c+1];
	if (l0 > r0) {
		*ll = 1;  *rr = 0;
		return -1;
	}
	r = csa->psi_pred(csa,r,l0,r0);
	l = csa->psi_succ(csa,l,l0,r0);
	if (l > r) {
		*ll = 1;  *rr = 0;
		return -1;
	}

	*ll = l;  *rr = r;
	return 0;
}

///////////////////////////////////////////
//  i64 csa_search(uchar *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri)
//    returns the length of the longest suffix of key using psi
//    [*li, *ri] is the range of key
///////////////////////////////////////////
i64 csa_search(uchar *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri)
{
	i64 c,h;
	i64 l,r,l0,r0,pl,pr;
	i64 len;

	//  pl= 1; pr = csa->n;
	pl= 0; pr = csa->n;
	if (keylen == 0) {
		*li = pl;  *ri = pr;
		return 0;
	}
	c = csa->CtoA[key[keylen-1]];
	r = csa->K[c+2]-1;  l = csa->K[c+1];
	len = 0;
	if (l > r) goto end;
	len++;
	for (h = keylen-2; h >= 0; h--) {
		pl = l;  pr = r;
		c = csa->CtoA[key[h]];
		r = csa->K[c+2]-1;  l = csa->K[c+1];
		if (l > r) goto end;

		l0 = l;  r0 = r;
		r = csa->psi_pred(csa,pr,l0,r0);
		l = csa->psi_succ(csa,pl,l0,r0);
		if (l > r) goto end;
		len++;
	}
end:
	if (len < keylen) {l = pl;  r = pr;}
	*li = l;  *ri = r;
	return len;
}

///////////////////////////////////////////
//  i64 csa_search_lf(uchar *key,i64 keylen,CSA *csa,rank_t *ll,rank_t *rr)
//    returns the length of the longest suffix of key using BW
//    [*ll, *rr] is the range of key
///////////////////////////////////////////
i64 csa_search_lf(uchar *key, i64 keylen, CSA *csa, rank_t *ll, rank_t *rr)
{
	int c;
	i64 l,r,len2,pl,pr;

	if (keylen == 0) {
		*ll = 0;  *rr = csa->n;
		return 0;
	}

	len2 = keylen;

	c = csa->CtoA[key[keylen-1]];
	l = csa->K[c+1];  r = csa->K[c+2]-1;
	//  printf("search\n");
	//  printf("[%ld,%ld]\n",l,r);
	while (--keylen > 0) {
		pl = l;  pr = r;
		c = csa->CtoA[key[keylen-1]];
		l = csa->K[c+1]-1 + csa->rankc(csa, l-1,c) + 1;
		r = csa->K[c+1]-1 + csa->rankc(csa, r,c);
		//    printf("[%ld,%ld]\n",l,r);
		if (l > r) goto end;
	}
end:
	if (keylen > 0) {l = pl;  r = pr;}
	*ll = l;  *rr = r;
	return len2-keylen;
}



///////////////////////////////////////////
//  i64 csa_search_r(i64 keylen,int c, CSA *csa,rank_t *li,rank_t *ri)
//    returns the length of the longest suffix of key:c
//      where key is the substring of length keylen represented by range [*li, *ri]
//    returned [*li, *ri] is the range of key:c
///////////////////////////////////////////
i64 csa_search_r(i64 keylen, int c, CSA *csa,rank_t *li,rank_t *ri)
{
	uchar *substr;
	i64 ret;

	substr = (uchar *) mymalloc(keylen+1);

	csa->substring(substr, csa, *li, keylen);
	substr[keylen] = c;

	ret = csa->search(substr, keylen+1, csa, li, ri);

	myfree(substr,keylen+1);

	return ret;
}


///////////////////////////////////////////
//  int csa_child_l(CSA *csa, rank_t l, rank_t r, uchar *head, rank_t *ll, rank_t *rr)
//    returns the number of characters preceding suffixes [l,r]
//  head stores the list of preceding characters
//  ll and rr store ranges for the preceding characters
//  lengths of head, ll, and rr must be at least csa->m (#distinct characters)
///////////////////////////////////////////
int csa_child_l(CSA *csa, rank_t l, rank_t r, uchar *head, rank_t *ll, rank_t *rr)
{
	int i,c,num;
	rank_t lll, rrr;
	i64 ret;

	num = 0;
	for (i=0; i<csa->m; i++) {
		c = csa->AtoC[i];
		lll = l;  rrr = r;
		ret = csa->searchsub(c, csa, &lll, &rrr);
		if (ret == 0) {
			head[num] = c;
			ll[num] = lll;
			rr[num] = rrr;
			num++;
		}
	}
	return num;
}

///////////////////////////////////////////
//  int csa_child_r(CSA *csa,i64 len,rank_t l,rank_t r,uchar *tail,rank_t *ll,rank_t *rr)
//    returns the number of characters following suffixes [l,r]
//  tail stores the list of following characters
//  ll and rr store ranges for the following characters
//  lengths of head, ll, and rr must be at least csa->m (#distinct characters)
///////////////////////////////////////////
int csa_child_r0(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr)
{
	int i,c,num;
	rank_t lll, rrr;
	uchar *substr;
	i64 ret;

	substr = (uchar *) mymalloc(len+1);
	csa->substring(substr, csa, l, len);

	num = 0;
	for (i=0; i<csa->m; i++) {
		c = csa->AtoC[i];
		lll = l;  rrr = r; // 不要?
		substr[len] = c;
		ret = csa->search(substr, len+1, csa, &lll, &rrr);
		if (ret == len+1) {
			tail[num] = c;
			ll[num] = lll;
			rr[num] = rrr;
			num++;
		}
	}

	myfree(substr,len+1);

	return num;
}

int csa_child_r_sub(CSA *csa, i64 len, rank_t l, rank_t r, int num, uchar *substr, uchar *tail, rank_t *ll, rank_t *rr)
{
	//int i,c;
	rank_t lll, rrr;
	i64 j, ret;
	int c1,c2;
	i64 l2 = 0, r2 = 0;

	//  printf("child_r_sub [%ld,%ld]\n",l,r);

	lll = l;  rrr = r;
	for (j=0; j<len; j++) {
		lll = csa->psi(csa, lll);
		rrr = csa->psi(csa, rrr);
	}
	c1 = csa->T(csa, lll);
	c2 = csa->T(csa, rrr);

	if (c1 == c2) {
		//    printf("num=%d c=%d [%ld,%ld]\n",num,c1,l,r);
		tail[num] = c1;
		ll[num] = l;
		rr[num] = r;
		num++;
	} else {
		if (c1 == -1) {
			l2 = l+1;
		} else {
			lll = l;  rrr = r;
			substr[len] = c1;
			ret = csa->search(substr, len+1, csa, &lll, &rrr);
			if (ret == len+1) {
				//        printf("num=%d c=%d [%ld,%ld]\n",num,c1,lll,rrr);
				tail[num] = c1;
				ll[num] = lll;
				rr[num] = rrr;
				num++;
				l2 = rrr+1;
			} else {
				printf("??? c1\n");
			}
		}
		if (c2 == -1) {
			r2 = r-1;
		} else {
			lll = l;  rrr = r;
			substr[len] = c2;
			ret = csa->search(substr, len+1, csa, &lll, &rrr);
			if (ret == len+1) {
				//        printf("num=%d c=%d [%ld,%ld]\n",num,c2,lll,rrr);
				tail[num] = c2;
				ll[num] = lll;
				rr[num] = rrr;
				num++;
				r2 = lll-1;
			} else {
				printf("??? c2\n");
			}
		}
		if (l2 <= r2) {
			num = csa_child_r_sub(csa, len, l2, r2, num, substr, tail, ll, rr);
		}
	}

	return num;
}

int csa_child_r(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr)
{
	int i,c,num;
	uchar *substr;

	rank_t *ltmp, *rtmp;
	uchar *tailtmp;
	int *hc;

	substr = (uchar *) mymalloc(len+1);
	csa->substring(substr, csa, l, len);

	tailtmp = (uchar *) mymalloc(SIGMA * sizeof(*tailtmp));
	ltmp = (rank_t *) mymalloc(SIGMA * sizeof(*ltmp));
	rtmp = (rank_t *) mymalloc(SIGMA * sizeof(*ltmp));
	hc = (int *) mymalloc(SIGMA * sizeof(*hc));

	num = csa_child_r_sub(csa, len, l, r, 0, substr, tailtmp, ltmp, rtmp);

	for (c=0; c<SIGMA; c++) hc[c] = -1;
	for (i=0; i<num; i++) hc[tailtmp[i]] = i;
	i = 0;
	for (c=0; c<SIGMA; c++) {
		if (hc[c] >= 0) {
			tail[i] = tailtmp[hc[c]];
			ll[i] = ltmp[hc[c]];
			rr[i] = rtmp[hc[c]];
			i++;
		}
	}

	free(hc);  free(rtmp);  free(ltmp);  free(tailtmp);

	myfree(substr,len+1);

	return num;
}

///////////////////////////////////////////
//  int csa_left_diverse(CSA *csa, rank_t l, rank_t r)
//    returns 1 if there exist two different preceding characters of a pattern (represented by rank [l,r]),
//    0 otherwise.
///////////////////////////////////////////
int csa_left_diverse(CSA *csa, rank_t l, rank_t r)
{
	int bw;
	rank_t ll, rr;

	bw = csa->BW(csa, l);
	ll = l;  rr = r;
	csa->searchsub(bw, csa, &ll, &rr);

	return (rr - ll != r - l);
}

///////////////////////////////////////////
//  int csa_right_diverse(CSA *csa, rank_t l, rank_t r, i64 length)
//    returns 1 if there exist two different following characters of a pattern (represented by rank [l,r] and length),
//    0 otherwise.
///////////////////////////////////////////
int csa_right_diverse(CSA *csa, rank_t l, rank_t r, i64 length)
{
	int c1,c2;
	i64 i;

	for (i=0; i<length; i++) {
		l = csa->psi(csa, l);
		r = csa->psi(csa, r);
	}

	c1 = csa->head(csa, l);
	c2 = csa->head(csa, r);

	return (c1 != c2);
}
