#include "ped_util.h"
/*
enum Condition get_condition_from_phenotype(int phenotype, ped_file_t *ped_file) {
	//if(ped_file->unaffected_id == -1)
        //ped_file->unaffected_id = kh_get(str, ped_file->phenotypes, "1");
    //if(ped_file->affected_id == -1)
        //ped_file->affected_id = kh_get(str, ped_file->phenotypes, "2");
    
    if(phenotype == -1){
        return UNKNOWN_CONDITION;
    }else if (phenotype == ped_file->unaffected_id) {//printf("ped_util.c : Phenotype %d es unaffected \n", phenotype);
        return UNAFFECTED;
    } else if (phenotype == ped_file->affected_id) {//printf("ped_util.c : Phenotype %d es affected\n", phenotype);
        return AFFECTED;
    } else {//printf("ped_util.c : Phenotype %d es desconocido\n", phenotype);
        return UNKNOWN_CONDITION;
    }
}*/

enum Condition get_condition_from_phenotype(char* phenotype, ped_file_t *ped_file) {
    if (!strcmp(phenotype, ped_file->unaffected)) {
        LOG_DEBUG_F("Phenotype %s is UNAFFECTED \n", phenotype);
        return UNAFFECTED;
    } else if (!strcmp(phenotype, ped_file->affected)) {
        LOG_DEBUG_F("Phenotype %s is AFFECTED\n", phenotype);
        return AFFECTED;
    } else {
        LOG_DEBUG_F("Phenotype %s is UNKNOWN\n", phenotype);
        return UNKNOWN_CONDITION;
    }
}


