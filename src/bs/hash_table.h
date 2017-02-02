#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define ZONE_CpG   1
#define ZONE_CHG   2
#define ZONE_CHH   3
#define ZONE_OTHER 0

typedef struct hash_node_s
{
    void* key;
    void* data;
    struct hash_node_s *next;
} hash_node_t;

typedef struct 
{
    hash_node_t *list;
    short count;
} hash_data_t;

typedef struct
{
    hash_data_t *table;
    int max_colitions;
    int size;
    int used;
    unsigned int (*hash_f)(void *k);
    int (*compare)(const void* a, const void* b); 
    void (*destroy_key)(void* a);
    void (*destroy_data)(void* a);
} hash_table_t;
