#include "ped_error.h"

char* get_ped_syntax_error_msg(enum PED_syntax_error error) {
    char *message = NULL;
    switch (error) {
        case COLUMN_FAMILY_ID_MISSING:
            message = "The column FAMILY_ID is missing";
            break;
        case COLUMN_INDIVIDUAL_ID_MISSING:
            message = "The column INDIVIDUAL_ID is missing";
            break;
        case COLUMN_FATHER_ID_MISSING:
            message = "The column FATHER_ID is missing";
            break;
        case COLUMN_MOTHER_ID_MISSING:
            message = "The column MOTHER_ID is missing";
            break;
        case COLUMN_SEX_MISSING:
            message = "The column SEX is missing";
            break;
        case COLUMN_PHENOTYPE_MISSING:
            message = "The column PHENOTYPE is missing";
            break;
        default:
            LOG_ERROR_F("Non-valid syntax error code (%d)\n", error);
    }
    return message;
}


char* get_ped_semantic_error_msg(enum PED_semantic_error error) {
    char *message = NULL;
    switch (error) {
        case ALREADY_EXISTING_FAMILY:
            message = "The family already exists";
            break;
        case CHILD_OF_VARIOUS_FAMILIES:
            message = "The individual is a child in more than one family";
            break;
        case NON_COHERENT_INDIVIDUAL_AMONG_FAMILIES:
            message = "The individual appears in several families and its data don't match";
            break;
        case FAMILY_WITH_MORE_THAN_TWO_FOUNDERS:
            message = "The family has more than two founders";
            break;
        case FATHER_WITH_FEMALE_SEX:
            message = "The individual is a father while having female sex";
            break;
        case MOTHER_WITH_MALE_SEX:
            message = "The individual is a mother while having male sex";
            break;
        case FATHER_APPEARS_MORE_THAN_ONCE:
            message = "The family already has a father set";
            break;
        case MOTHER_APPEARS_MORE_THAN_ONCE:
            message = "The family already has a mother set";
            break;
        default:
            LOG_ERROR_F("Non-valid semantic error code (%d)\n", error);
    }
    return message;
}
