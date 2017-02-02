#include "region_table.h"

region_table_t *create_table(const char *url, const char *species, const char *version)
{
    int num_chromosomes;
    region_table_t *table = (region_table_t*) malloc (sizeof(region_table_t));
	
    table->ordering = get_chromosome_order(url, species, version, &num_chromosomes);
    table->max_chromosomes = num_chromosomes;

    table->storage = cp_hashtable_create_by_option( COLLECTION_MODE_NOSYNC | COLLECTION_MODE_DEEP,
                                                    num_chromosomes * 2,
                                                    cp_hash_istring,			// Hash function
                                                    (cp_compare_fn) strcasecmp,		// Key comparison function
                                                    NULL,					// Key copy function
                                                    (cp_destructor_fn) free_chromosome,	// Key destructor function
                                                    NULL,					// Value copy function
                                                    (cp_destructor_fn) free_chromosome_tree	// Value destructor function
    );

    return table;
}

void free_table(region_table_t* regions) {
    // Free ordering array
    char **ordering = regions->ordering;
    for (int i = 0; i < regions->max_chromosomes; i++) {
        free(ordering[i]);
    }
    free(ordering);

    // Free hashtable
    cp_hashtable *table = regions->storage;
    cp_hashtable_destroy(table);
    free(regions);
}


/* ******************************
 *  		Regions		*
 * ******************************/

int insert_region(region_t *region, region_table_t *table)
{
	cp_avltree *chr = NULL;
	if (!cp_hashtable_contains(table->storage, region->chromosome))
	{
        // Check if the chromosome is valid for the species
        int valid = 0;
        for (int i = 0; i < table->max_chromosomes && !valid; i++) {
            if (!strcasecmp(region->chromosome, table->ordering[i])) {
                valid = 1;
            }
        }
        if (!valid) {
            LOG_WARN_F("Chromosome %s does not match the species to analyze\n", region->chromosome);
            return -1;
        }
        
		// Insert new chromosome
		if ((chr = insert_chromosome(region->chromosome, table)) == NULL)
		{
			return 1;
		}
	}

	// Insert region in chromosome tree
	if (chr == NULL) 
	{
		chr = get_chromosome(region->chromosome, table);
	}
	if (cp_avltree_insert(chr, &region->start_position, region) == NULL)
	{
		return 2;
	}
	return 0;
}

int contains_region(region_t *region, region_table_t *table)
{
	cp_avltree *chr = NULL;
	if (region == NULL || table == NULL || 
		(chr  = get_chromosome(region->chromosome, table)) == NULL)
// 		!cp_hashtable_contains(table->storage, region->chromosome))
	{
		return 0;
	}
	
	return cp_avltree_contains(chr, &(region->start_position));
}

int find_region(region_t *region, region_table_t *table)
{
	cp_avltree *chr = NULL;
	if (region == NULL || table == NULL || 
		(chr  = get_chromosome(region->chromosome, table)) == NULL)
// 		!cp_hashtable_contains(table->storage, region->chromosome))
	{
		return 0;
	}
	
	// Set new comparison function
	cp_compare_fn aux_fn = (cp_compare_fn) chr->cmp;
	chr->cmp = (cp_compare_fn) region_contains_other;
	// Check whether it contains the region
	int result = cp_avltree_contains(chr, region);
	// Restore comparison function
	chr->cmp = aux_fn;

	return result;
}

region_t *remove_region(region_t *region, region_table_t *table)
{
	if (!cp_hashtable_contains(table->storage, region->chromosome))
	{
		return NULL;
	}
	
	region_t *removed = NULL;
	cp_avltree *chr = get_chromosome(region->chromosome, table);
	/*
	 * - If the tree node contains >1 element: remove all elements in that exact range. 
	 *   If all of them are erased, remove the node.
	 * - If the tree node contains 1 element: remove the node
	 */
	cp_vector *region_vector = cp_avltree_get(chr, &region->start_position);
	int regions_in_node = cp_vector_size(region_vector);
	
	if (regions_in_node > 1)
	{
		int i = 0;
		do
		{
			// Remove all elements in that exact range. 
			if (compare_position_ranges(region, cp_vector_element_at(region_vector, i)) == 0)
			{
				removed = cp_vector_remove_element_at(region_vector, i);
				regions_in_node--;
			} else
			{
				i++;
			}
		} while (i < regions_in_node);
		
		// If all elements are erased, remove the node
		regions_in_node = cp_vector_size(region_vector);
		if (regions_in_node == 0)
		{
			removed = cp_avltree_delete(chr, &region->start_position);
		}
	} else
	{
		removed = cp_avltree_delete(chr, &region->start_position);
	}
	
	return removed;
}

void free_region(region_t *region)
{
	// Doesn't need to be dynamically allocated memory
// 	free(region->chromosome);
	free(region);
}



/* ******************************
 *  Chromosomes (region trees)	*
 * ******************************/

void *insert_chromosome(const char *key, region_table_t *table)
{
	cp_avltree *chr_tree = get_chromosome(key, table);
	if (chr_tree != NULL)
	{
		return chr_tree;
	}
	chr_tree = cp_avltree_create_by_option( COLLECTION_MODE_MULTIPLE_VALUES | COLLECTION_MODE_DEEP,
// 						(cp_compare_fn) compare_regions,	// Key comparison function
						(cp_compare_fn) compare_positions,
						NULL,					// Key copy function
						NULL,					// Key destructor function
						NULL,					// Value copy function
						(cp_destructor_fn) free_region		// Value destructor function
						);
	return cp_hashtable_put(table->storage, (void*) key, chr_tree);
}

cp_avltree *get_chromosome(const char *key, region_table_t *table)
{
	return (cp_avltree *) cp_hashtable_get_by_option(table->storage, (void*) key, COLLECTION_MODE_NOSYNC);
}

cp_avltree *remove_chromosome(const char *key, region_table_t *table)
{
	return cp_hashtable_remove(table->storage, (void*) key);
}

int count_regions_in_chromosome(const char *key, region_table_t *table) {
    cp_avltree *chr_tree = get_chromosome(key, table);
    if (chr_tree == NULL) {
        return 0;
    }
    return cp_avltree_count(chr_tree);
}

void free_chromosome(const char *key)//, region_table_t *table)
{
	// As in free_region, the key doesn't need to be dynamically allocated memory
// 	free((void*) key);
}

void free_chromosome_tree(cp_avltree *chromosome_tree)
{
	cp_avltree_destroy(chromosome_tree);
}
