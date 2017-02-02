#ifndef PED_ERROR_H
#define PED_ERROR_H

#include <commons/log.h>

enum PED_syntax_error {
    COLUMN_FAMILY_ID_MISSING = 100,
    COLUMN_INDIVIDUAL_ID_MISSING,
    COLUMN_FATHER_ID_MISSING,
    COLUMN_MOTHER_ID_MISSING,
    COLUMN_SEX_MISSING,
    COLUMN_PHENOTYPE_MISSING
};


enum PED_semantic_error {
    ALREADY_EXISTING_FAMILY = 200,
    CHILD_OF_VARIOUS_FAMILIES,
    NON_COHERENT_INDIVIDUAL_AMONG_FAMILIES,
    FAMILY_WITH_MORE_THAN_TWO_FOUNDERS,
    FATHER_WITH_FEMALE_SEX,
    MOTHER_WITH_MALE_SEX,
    FATHER_APPEARS_MORE_THAN_ONCE,
    MOTHER_APPEARS_MORE_THAN_ONCE,
};


char *get_ped_syntax_error_msg(enum PED_syntax_error error);

char *get_ped_semantic_error_msg(enum PED_semantic_error error);

#endif