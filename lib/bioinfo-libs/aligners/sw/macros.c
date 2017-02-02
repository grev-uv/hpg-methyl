#include "macros.h"

//------------------------------------------------------------------------
// utility functions
//------------------------------------------------------------------------

void str_trim(char *line) {
  int len = strlen(line) - 1 ;
  if (line[len] == '\n' || line[len] == '\r') {
    line[len] = '\0';
  }
}

//------------------------------------------------------------------------

double sw_tic() {
  struct timeval t1;
  gettimeofday(&t1, NULL);
  return ((double) t1.tv_sec*1000000 + (double) t1.tv_usec) / 1000000.0f;
}

//------------------------------------------------------------------------

double sw_toc(double t1) {
  struct timeval t2;

  gettimeofday(&t2, NULL);
  return (((double) t2.tv_sec*1000000 + (double) t2.tv_usec) / 1000000.0f) - t1;
}


//------------------------------------------------------------------------

void simd_traceback(int depth, int num_seqs, 
		    char **q, int *q_len, int max_q_len,
		    char **r, int *r_len, int max_r_len,
		    float gap_open, float gap_extend,
		    float *H, int *C,
		    float *max_score,
		    char **q_alig, int *q_start,
		    char **r_alig, int *r_start, 
		    int *len_alig,
		    char *q_aux, char *r_aux) {
  char *qq, *qq_alig;
  char *rr, *rr_alig;
  
  int jj, kk;

  int ix, iy, gapcnt, index;
  double bimble;
  double errbounds = (double) 0.01;

  int qq_len, rr_len;
  int max_len = 0, len = 0;
  char c;
  float score, diagonal, e, f;

  for (int i = 0; i < num_seqs ; i++) {
    qq = q[i];
    rr = r[i];
    qq_len = q_len[i];
    rr_len = r_len[i];

    simd_find_position(depth, i, qq, qq_len, rr, rr_len, H, max_q_len, max_r_len, max_score[i], &kk, &jj);
    len = 0;

    while (1) {
      index = (jj * max_q_len * depth) + (kk * depth) + i;
      score = H[index];
      
      if (score == 0.0f) {
        break;
      }

      c = C[index];

      if (c == 5 || c == 7) { // diagonal
        q_aux[len] = qq[kk];
        r_aux[len] = rr[jj];

        len++;
        jj--;
        kk--;

        if (jj < 0) break;
        if (kk < 0) break;

        continue;
      } else if (c == 1 || c == 3) { // left
        gapcnt = 0;
        ix     = kk - 1;

        while (1) {
          bimble = H[(jj * max_q_len * depth) + (ix * depth) + i] - gap_open - (gapcnt * gap_extend);

          if (!ix || fabs((double)score - (double)bimble) < errbounds) {
            break;
          }

          --ix;

          if (ix < 0) { 
            printf("Error reconstructing alignment (walking left)\n");
            exit(-1);
          }

          ++gapcnt;
        }

        if(bimble <= 0.0)
        break;

        for(int ic = 0; ic <= gapcnt; ++ic) {
        q_aux[len] = qq[kk--]; 
        r_aux[len] = '-';
        len++;
        }

        continue;
      } else if (c == 2 || c == 6) { // down
        gapcnt = 0;
        iy     = jj - 1;

        while (1) {
          bimble = H[(iy * max_q_len * depth) + (kk * depth) + i] - gap_open - (gapcnt * gap_extend);

          if (!iy || fabs((double)score - (double)bimble) < errbounds) {
            break;
          }

          --iy; 

          if (iy < 0) { 
            printf("Error reconstructing alignment (walking down)\n");
            exit(-1);
          }

          ++gapcnt;
        }

        if(bimble <= 0.0) {
          break;
        }

        for(int ic = 0; ic <= gapcnt; ++ic) {
          q_aux[len] = '-';
          r_aux[len] = rr[jj--];
          len++;
        }

        continue;
      } else {
        break;
      }
    }


    q_start[i] = kk + 1;
    r_start[i] = jj + 1;
    
    q_aux[len] = 0;
    r_aux[len] = 0;

    qq_alig = (char *) calloc(len + 1, sizeof(char));
    rr_alig = (char *) calloc(len + 1, sizeof(char));
    jj = len - 1;

    for (int i = 0; i < len; i++) {
      qq_alig[i] = q_aux[jj];
      rr_alig[i] = r_aux[jj];
      jj--;
    }

    qq_alig[len] = 0;
    rr_alig[len] = 0;

    q_alig[i] = qq_alig;
    r_alig[i] = rr_alig;

    len_alig[i] = len;
  }
}

//-------------------------------------------------------------------------

void simd_find_position(int depth, int index, char *q, int q_len, char *r, int r_len,
			float *H, int cols, int rows, float score, 
			int *q_pos, int *r_pos) {  
  *r_pos = 0;
  *q_pos = 0;
  int i, ii = 0;

  for (int j = 0; j < r_len; j++) {
    i = j * cols * depth;
    
    for (int k = 0; k < q_len; k++) {
      ii = i + (k * depth) + index;
      
      if (H[ii] == score) {
        *r_pos = j;
        *q_pos = k;
        return;
      }
    }
  }
}

//------------------------------------------------------------------------

void save_float_matrix(float *matrix, int cols, int rows, char *q, int q_len, 
		       char *r, int r_len, int index, int depth, char *filename) {
  
  FILE *f = fopen(filename, "w");

  if (f == NULL) {
    return;
  }

  fprintf(f, "\t");

  for (int j = 0; j < q_len; j++) {
    fprintf(f, "%c\t", q[j]);
  }
  
  fprintf(f, "\n");

  for (int i = 0; i < r_len; i++) {
    fprintf(f, "%c\t", r[i]);

    for (int j = 0; j < q_len; j++) {
      fprintf(f, "%0.2f\t", matrix[(i * cols * depth) + (j * depth) + index]);
    }

    fprintf(f, "\n");
  }
  
  fclose(f);
}

//------------------------------------------------------------------------

void save_int_matrix(int *matrix, char *q, int q_len, char *r, int r_len, 
		     int index, int depth, char *filename) {
  FILE *f = fopen(filename, "w");

  if (f == NULL) {
    return;
  }

  fprintf(f, "\t");

  for (int k = 0; k < q_len; k++) {
    fprintf(f, "%c\t", q[k]);
  }

  fprintf(f, "\n");

  int i, ii;

  for (int j = 0; j < r_len; j++) {
    fprintf(f, "%c\t", r[j]);
    i = j * q_len * depth;

    for (int k = 0; k < q_len; k++) {
      ii = i + (k * depth) + index;
      fprintf(f, "%i\t", matrix[ii]);
    }

    fprintf(f, "\n");
  }

  fclose(f);
}
