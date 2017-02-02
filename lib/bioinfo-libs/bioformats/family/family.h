#ifndef FAMILY_H
#define FAMILY_H

#include <assert.h>
#include <stdlib.h>


#include <commons/log.h>
#include <commons/string_utils.h>
#include <containers/linked_list.h>


enum Sex { MALE, FEMALE, UNKNOWN_SEX };

/**
 * Missing: Created because one of his child was read before, but this individual has no information available.
 * Affected: This individual is affected by a disease.
 * Unaffected: This individual is not affected by a disease.
 * Unknown: The condition of the individual is unknown (phenotype value is not valid).
 */
enum Condition { MISSING_CONDITION, AFFECTED, UNAFFECTED, UNKNOWN_CONDITION };

/**
 * Entry in the PED document body, representing an individual and member of a family.
 */
typedef struct individual {
    char *id;   /**< Unique ID of the individual **/
    float phenotype;    /**< Numerical descriptor for the affection of the individual */
    enum Sex sex;   /**< Sex of the individual */
    enum Condition condition;  /**< Whether the individual is affected by a disease or not, or its data are missing */
    struct individual *father;  /**< Father of the individual (NULL in case he's parent in a family) */
    struct individual *mother;  /**< Mother of the individual (NULL in case he's parent in a family) */
    struct family *family;  /**< Family of the individual **/
} individual_t;

/**
 * Family described in a PED document
 */
typedef struct family {
    char *id;               /**< Unique ID of the family **/
    individual_t *father;   /**< Man in the root of the genealogical tree */
    individual_t *mother;   /**< Woman in the root of the genealogical tree */
    linked_list_t *children;      /**< Children of the main roots in the genealogical tree */
    linked_list_t *unknown;       /**< Unclassified samples because they have no parents and no sex */
} family_t;


/*
 * Individual management functions
 */

/**
 * Creates a new individual with the characteristics provided as arguments.
 * 
 * @param id unique ID of the individual
 * @param phenotype numerical descriptor for the affection of the individual
 * @param sex sex of the individual
 * @param condition whether the individual is un/affected or missing
 * @param father father of the individual, if apply
 * @param mother mother of the individual, if apply
 * @param family family of the individual
 * @return The newly created individual
 */
individual_t *individual_new(char *id, float phenotype, enum Sex sex, enum Condition condition, individual_t *father, individual_t *mother, family_t *family);

/**
 * Fills member of an already existing individual with the characteristics provided as arguments.
 * 
 * @param id unique ID of the individual
 * @param phenotype numerical descriptor for the affection of the individual
 * @param sex sex of the individual
 * @param father father of the individual, if apply
 * @param mother mother of the individual, if apply
 * @param family family of the individual
 * @param individual individual whose members are being filled
 */
void individual_init(char *id, float phenotype, enum Sex sex, enum Condition condition, individual_t *father, individual_t *mother, family_t *family, individual_t *individual);

/**
 * Free memory associated to an individual.
 * 
 * @param individual individual to be freed
 */
void individual_free(individual_t *individual);


/*
 * Family management functions
 */

/**
 * Creates a new family with the characteristics provided as arguments.
 * 
 * @param id unique ID of the family
 * @return The newly created family
 * 
 */
family_t *family_new(char *id);

/**
 * Sets an individual as parent of a family (its position as father or mother is automatically 
 * assigned based on its sex). An individual can only be set as parent if his sex is known and 
 * he is not part of the family yet.
 * 
 * @param parent the individual to set as parent
 * @param family the family where the individual could belong to
 * @return 0 if the individual was succcessfully set, 1-2 if one of the arguments is NULL, 
 * 3 if the parent sex in unknown, or 4 if the family already contains the individual
 */
int family_set_parent(individual_t *parent, family_t *family);

/**
 * Sets an individual as child of a family. An individual can only be set as child if he is 
 * not part of the family yet.
 * 
 * @param child the individual to set as child
 * @param family the family where the individual could belong to
 * @return 0 if the individual was succcessfully set, 1-2 if one of the arguments is NULL, 
 * 3 if the family already contains the individual
 */
int family_add_child(individual_t *child, family_t *family);


int family_add_unknown(individual_t *individual, family_t *family);


/**
 * Free memory associated to a family and its individuals.
 * 
 * @param family the family to be freed
 */
void family_free(family_t *family);

/**
 * Checks whether an individual is part of a family.
 * 
 * @param individual the individual to check
 * @param family the family that could contain the individual
 */
individual_t *family_contains_individual(individual_t *individual, family_t *family);

#endif
