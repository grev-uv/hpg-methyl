#ifndef REGION_H
#define REGION_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include <curl/curl.h>

#include <commons/http_utils.h>
#include <commons/log.h>
#include <commons/string_utils.h>

typedef struct region {
    size_t start_position;
    size_t end_position;
    char *chromosome;
} region_t;

typedef struct {
    char *data;
    size_t length;
} chromosome_ws_response;

region_t *region_new(char *chromosome, uint32_t start_position, uint32_t end_position);

void region_free(region_t *region);

/**
 * Get an ordered list of chromosomes which will be used as reference in 
 * regions' comparison functions.
 * 
 * If a file is provided its order is used, otherwise the default criteria
 * is inserting the numerical chromosomes at the head of the list, then
 * use lexicographical order for textual names.
 * 
 * @param species
 * 	(Optional) File the list of ordered chromosomes is read from
 * @param num_chromosomes
 * 	(Output) Number of chromosomes present in the list
 * 
 * @return
 * 	List containing the names of the chromosomes in order
 */
char **get_chromosome_order(const char *host_url, const char *species, const char *version, int *num_chromosomes);

static char *compose_chromosomes_ws_request(const char *host_url, const char *species, const char *version);

static size_t write_chromosomes_ws_results(char *contents, size_t size, size_t nmemb, chromosome_ws_response *userdata);

/**
 * Compare two regions, considering their chromosome and start position.
 * 
 * @param region_1
 * 	First region to compare
 * @param region_2
 * 	Second region to compare
 * @param chromosome_ordering
 * 	Chromosome ordering criteria
 * @param num_chromosomes
 * 	Number of chromosomes contained in 'chromosome_ordering'
 * 
 * @return
 * 	Less, equals or greater than zero if region_1 is placed before, in the same position or after region_2.
 * 	INT_MIN if any of the regions is NULL.
 */
int compare_regions(void *region_1, void *region_2, char **chromosome_ordering, int num_chromosomes);

int compare_chromosomes(char *chromosome_1, char *chromosome_2, char **chromosome_ordering, int num_chromosomes);

/**
 * Compare two regions' positions (could be start or end).
 * 
 * @param position_1
 * 	Position of the first region to compare
 * @param position_2
 * 	Position of the second region to compare
 * 
 * @return
 * 	Less, equal or greater than zero if position_1 is less than, equal or greater than position_2.
 */
int compare_positions(uint32_t position_1, uint32_t position_2);

/**
 * Compare two regions' position ranges.
 * 
 * @param region_1
 * 	First region to compare
 * @param region_2
 * 	Second region to compare
 * 
 * @return
 * 	Less, equal or greater than zero if position_1 is less than, equal or greater than position_2.
 */
int compare_position_ranges(region_t *region_1, region_t *region_2);


/**
 * Check whether a region called 'content' is contained inside other called 'container'. One region 
 * contains another if both start and end positions of the former cover, at least, the whole range of the last.
 * 
 * @param container
 * 	First region to compare
 * @param region_2
 * 	Second region to compare
 * 
 * @return
 * 	Less, equal or greater than zero if the content is positioned before, inside or after the container.
 */
int region_contains_other(region_t *container, region_t *content);

#endif