#include "rna_map_region.h"


void rna_map_region_free(void* region)
{
	rna_map_region_t *reg = (rna_map_region_t*) region;
	free(reg->chromosome);
	free(reg);
}

int rna_map_region_equal(rna_map_region_t *region1,rna_map_region_t* region2){
	int result=0;
	if(region1->strand != region2->strand){
		result = ((region1->start_position == region2->end_position) || (region2->start_position == region1->end_position))
			&& (strcmp(region1->chromosome,region2->chromosome) == 0);
		// 	 region1->end_position == region2->start_position;

	}
	else{
		result = ( region1->start_position == region2->start_position
				)						&& (strcmp(region1->chromosome,region2->chromosome) == 0 );
		// region1->end_position == region2->end_position;
		//strcaseicmp(region1->chromosome,region2->chromosome)==0 
	}

	return result;
}

void rna_print_region(rna_map_region_t* region){
	printf("Ch: %s - %d - %d - %c",region->chromosome, region->start_position,region->end_position,region->strand == 0 ? '0':'1' );
}