#ifndef REGION_TABLE_UTILS_H
#define REGION_TABLE_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <limits.h>
#include <omp.h>

#include <commons/log.h>

#include <bioformats/bed/bed_file.h>
#include <bioformats/bed/bed_file_structure.h>
#include <bioformats/gff/gff_file.h>
#include <bioformats/gff/gff_file_structure.h>
#include <bioformats/features/region/region.h>
#include <bioformats/features/region/region_table.h>

region_table_t *parse_regions(char *input_regions, int as_positions, const char *url, const char *species, const char *version);

region_table_t *parse_regions_from_gff_file(char *filename, const char *url, const char *species, const char *version);

region_table_t *parse_regions_from_bed_file(char *filename, const char *url, const char *species, const char *version);

int region_table_parse_from_gff_file(char *filename, region_table_t *regions_table);
int region_table_parse_from_string(char *input_regions, region_table_t *regions_table);

#ifdef __cplusplus
}
#endif

#endif
