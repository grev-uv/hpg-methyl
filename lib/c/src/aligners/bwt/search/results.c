#include "results.h"

bool write_results(results_list *r_list, intmax_t *k, intmax_t *l, exome* ex, bwt_index *backward, bwt_index *forward, char *mapping, uintmax_t nW, int type, FILE *fp, bwt_config_t *bwt_config) {

	result *r;

	bool found = false;

	char search[MAXLINE+1];

  search[0] = '\0';
  strncat(search, mapping, nW);

	bool repeated;
  uintmax_t kl_count = 0;

  for (uintmax_t i=0;i<r_list->num_results; i++) {

		r = &r_list->list[i];

		repeated = false;
		uintmax_t kl;
		for (kl=0; kl < kl_count; kl++) {

			if(r->k == k[kl] && r->l == l[kl]) {
				repeated = true;
				break;
			}

		}

		if (!repeated) {
			k[kl] = r->k;
			l[kl] = r->l;
			kl_count++;
			manage_single_result(r, ex, backward, forward, search, type, fp, r_list->read_index, &found, bwt_config);
		}

	}

	return found;

}
