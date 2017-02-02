#include "region.h"


region_t *region_new(char *chromosome, uint32_t start_position, uint32_t end_position) {
    assert(chromosome);
    region_t *region = malloc(sizeof(region));
    region->chromosome = chromosome;
    region->start_position = start_position;
    region->end_position = end_position;
    return region;
}

void region_free(region_t *region) {
    assert(region);
    free(region->chromosome);
    free(region);
}



char **get_chromosome_order(const char *host_url, const char *species, const char *version, int *num_chromosomes)
{
    int ret_code = init_http_environment(0);
    if (ret_code != 0) {
        return NULL;
    }
    
    char **ordering = NULL;
    
    // Default species: hsa
    if (species == NULL || version == NULL) {
        *num_chromosomes = 25;
        ordering = (char**) malloc ((*num_chromosomes) * sizeof(char*));
        for (int i = 0; i < 22; i++)
        {
            ordering[i] = (char*) calloc (3, sizeof(char));
            sprintf(ordering[i], "%d", i+1);
        }
        ordering[22] = (char*) calloc (2, sizeof(char));
        ordering[23] = (char*) calloc (2, sizeof(char));
        ordering[24] = (char*) calloc (3, sizeof(char));
        strcat(ordering[22], "X");
        strcat(ordering[23], "Y");
        strcat(ordering[24], "MT");
        
        return ordering;
    }
    
    CURL *curl;
    CURLcode res;
    char *url;

    curl = curl_easy_init();
    if(curl) {
        url = compose_chromosomes_ws_request(host_url, species, version);
        if (!url) { return NULL; }

        chromosome_ws_response s;
        s.length = 0;
        s.data = calloc(1, sizeof(char));

        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_chromosomes_ws_results);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &s);
        res = curl_easy_perform(curl);

        LOG_DEBUG_F("Species chromosomes = { %s }\n", s.data);
        ordering = split(s.data, ",", num_chromosomes);

        free(s.data);

        /* always cleanup */
        curl_easy_cleanup(curl);
        
        free(url);
    }
    
    return ordering;
}

static char *compose_chromosomes_ws_request(const char *host_url, const char *species, const char *version) {
    if (host_url == NULL || version == NULL || species == NULL) {
        return NULL;
    }
    
    // URL Constants
    const char *ws_root_url = "cellbase/rest/";
    const char *ws_name_url = "chromosomes";
    
    // Length of URL parts
    const int host_url_len = strlen(host_url);
    const int ws_root_len = strlen(ws_root_url);
    const int version_len = strlen(version);
    const int species_len = strlen(species);
    const int ws_name_len = strlen(ws_name_url);
    const int result_len = host_url_len + ws_root_len + version_len + species_len + ws_name_len + 4;
    
    char *result_url = (char*) calloc (result_len, sizeof(char));
    
    // Host URL
    strncat(result_url, host_url, host_url_len);
    if (result_url[host_url_len - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Root of the web service
    strncat(result_url, ws_root_url, ws_root_len);
    
    // Version
    strncat(result_url, version, version_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Species
    strncat(result_url, species, species_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Name of the web service
    strncat(result_url, ws_name_url, ws_name_len);
    
    return result_url;
}

static size_t write_chromosomes_ws_results(char *contents, size_t size, size_t nmemb, chromosome_ws_response *s) {
    size_t new_len = s->length + size*nmemb;
    s->data = realloc(s->data, new_len+1);
    if (s->data == NULL) {
        LOG_FATAL("Can't allocate enough memory for getting chromosomes");
    }
    memcpy(s->data+s->length, contents, size*nmemb);
    s->data[new_len] = '\0';
    s->length = new_len;

    return size*nmemb;
}



int compare_regions(void *region_1, void *region_2, char **chromosome_ordering, int num_chromosomes)
{
	if (region_1 == NULL || region_2 == NULL)
	{
		return INT_MIN;
	}
	
	region_t *reg_1 = (region_t *) region_1;
	region_t *reg_2 = (region_t *) region_2;
	
	// TODO This could be avoided while inserting, because regions are classified by chromosome
	int result = compare_chromosomes(reg_1->chromosome, reg_2->chromosome, chromosome_ordering, num_chromosomes);
	if (result != 0)
	{
		return result;
	} else
	{
// 		return compare_position_ranges(reg_1, reg_2);
		return compare_positions(reg_1->start_position, reg_2->start_position);
	}
}

int compare_chromosomes(char *chromosome_1, char *chromosome_2, char **chromosome_ordering, int num_chromosomes)
{
    assert(chromosome_1);
    assert(chromosome_2);
//     printf("chr1 = %s\t", chromosome_1);
//     printf("chr 2 = %s\t", chromosome_2);
//     printf("num chr = %d\n", num_chromosomes);
	int chr_1_found = 0, chr_2_found = 0;
	for (int i = 0; i < num_chromosomes; i++)
	{
        assert(chromosome_ordering[i]);
//         printf("* 0\n");
		if (strcasecmp(chromosome_ordering[i], chromosome_1) == 0)
		{
//             printf("* 1\n");
			if (strcasecmp(chromosome_ordering[i], chromosome_2) == 0)
			{
//                 printf("* 2\n");
				return 0;
			}
			return -1;
		}
		else if (strcasecmp(chromosome_ordering[i], chromosome_2) == 0)
		{
//             printf("* 4\n");
			return 1;
		}
	}
	return 0;
}

int compare_positions(uint32_t position_1, uint32_t position_2)
{
	return position_1 - position_2;
}

int compare_position_ranges(region_t *region_1, region_t *region_2)
{
	int result = region_1->start_position - region_2->start_position;
	if (result == 0)
	{
		result = region_1->end_position - region_2->end_position;
	}
	return result;
}

int region_contains_other(region_t *container, region_t *content)
{
// 	printf("container = %s:%d:%d\t", container->chromosome, container->start_position, container->end_position);
// 	printf("content = %s:%d:%d\n", content->chromosome, content->start_position, content->end_position);
// 	printf("start ok = %d, end ok = %d\n", container->start_position <= content->start_position, container->end_position >= content->end_position);
	
	int result = strcasecmp(container->chromosome, content->chromosome);
	if (result != 0) { return result; }
	
	if (container->start_position > content->start_position)
	{
		return 1;
	}
	
	if (container->end_position < content->end_position)
	{
		return -1;
	}
	
	return 0;
}
