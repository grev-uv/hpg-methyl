#include "dna_map_region.h"


void dna_map_region_free(void* region)
{
	dna_map_region_t *reg = (dna_map_region_t*) region;
	free(reg->chromosome);
	free(reg);
}

int dna_map_region_equal_hard(dna_map_region_t *region1,dna_map_region_t* region2){
	int result=0;
	if(region1->strand != region2->strand){
		result = (region1->start_position == region2->end_position) && 
				 (region2->start_position == region1->end_position) &&
				 (strcmp(region1->chromosome,region2->chromosome) == 0);

	}
	else{
		result = (region1->start_position == region2->start_position) && 
		         (region1->end_position == region2->end_position) &&
			     (strcmp(region1->chromosome,region2->chromosome) == 0 );
				  
	}

	return result;
}

int dna_map_region_equal_soft(dna_map_region_t *region1,dna_map_region_t* region2){
	int result=0;
	if(region1->strand != region2->strand){
		result = ((region1->start_position == region2->end_position) || 
				 (region2->start_position == region1->end_position)) &&
				 (strcmp(region1->chromosome,region2->chromosome) == 0);

	}
	else{
		result = (region1->start_position == region2->start_position)  &&
			     (strcmp(region1->chromosome,region2->chromosome) == 0 );
				  
	}

	return result;
}


void dna_print_region(dna_map_region_t* region){
	printf("Ch: %s - %d - %d - %c",region->chromosome, region->start_position,region->end_position,region->strand == 0 ? '0':'1' );
}

void dna_fprint_region(FILE *f,dna_map_region_t* region){
	fprintf(f,"Ch: %s - %d - %d - %c\n",region->chromosome, region->start_position,region->end_position,region->strand == 0 ? '0':'1' );
}
