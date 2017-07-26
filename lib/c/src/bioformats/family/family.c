#include "family.h"


/*
 * Individual management functions
 */

individual_t *individual_new(char *id, float variable, enum Sex sex, enum Condition condition, 
                             individual_t *father, individual_t *mother, family_t *family) {
    individual_t *individual = (individual_t*) malloc (sizeof(individual_t));
    individual_init(id, variable, sex, condition, father, mother, family, individual);
    return individual;
}

individual_t *individual_new_ids_only(char *id, float variable, enum Sex sex, enum Condition condition, 
                             char *father_id, char *mother_id, family_t *family) {
    individual_t *individual = (individual_t*) malloc (sizeof(individual_t));
    individual_init_ids_only(id, variable, sex, condition, father_id, mother_id, family, individual);
    return individual;
}


void individual_init(char *id, float variable, enum Sex sex, enum Condition condition, 
                    individual_t *father, individual_t *mother, family_t *family, individual_t *individual) {
    assert(individual);
    individual->id = id;
    individual->variable = variable;
    individual->condition = condition;
    individual->sex = sex;
    individual->father_id = father->id;
    individual->mother_id = mother->id;
    individual->father = father;
    individual->mother = mother;
    individual->children = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    individual->family = family;
}

void individual_init_ids_only(char *id, float variable, enum Sex sex, enum Condition condition, 
                              char *father_id, char *mother_id, family_t *family, individual_t *individual) {
    assert(individual);
    individual->id = id;
    individual->variable = variable;
    individual->condition = condition;
    individual->sex = sex;
    individual->father_id = father_id;
    individual->mother_id = mother_id;
    individual->father = NULL;
    individual->mother = NULL;
    individual->children = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    individual->family = family;
}

void individual_free(individual_t *individual) {
    assert(individual);
    free(individual->id);
    free(individual->father_id);
    free(individual->mother_id);
    linked_list_free(individual->children, NULL);
    free(individual);
}

int individual_add_child(individual_t *child, individual_t *individual) {
    return !linked_list_insert(child, individual->children);
}

int individual_compare(individual_t *a, individual_t *b) {
    int result;
    if (!a || !b) {
        result = 1;
    } else if (a->family == NULL && b->family != NULL) {
        result = -1;
    } else if (a->family != NULL && b->family == NULL) {
        result = 1;
    } else {
        result = strcasecmp(a->id, b->id);
        LOG_DEBUG_F("1) %s:%s = %s:%s ? %d\n", a->family->id, a->id, b->family->id, b->id, !result);
        result |= strcasecmp(a->family->id, b->family->id);
        LOG_DEBUG_F("2) %s:%s = %s:%s ? %d\n", a->family->id, a->id, b->family->id, b->id, !result);
    }
    return result;
}


/*
 * Family management functions
 */

family_t *family_new(char *id) {
    family_t *family = (family_t*) calloc (1, sizeof(family_t));
    family->id = id;
    family->members = kh_init(family_members);
    family->founders = kh_init(family_members);
    family->unknown = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    
    return family;
}

static int family_add_founder(individual_t *individual, family_t *family) {
    assert(individual);
    assert(family);
   
    int ret = 0;
    khiter_t iter = kh_put(family_members, family->founders, strdup(individual->id), &ret);
    if (ret) {
        kh_value(family->founders, iter) = individual;
        return 0;
    }
    return 1;
}

int family_add_member(individual_t *individual, family_t *family) {
    assert(individual);
    assert(family);
    
    if (family_contains_individual(individual, family)) {
        // TODO return correct error code
        return -1;
    }
    
    int ret = 0;
    // Combine the exit codes of adding the individual to the members...
    khiter_t iter = kh_put(family_members, family->members, strdup(individual->id), &ret);
    if (ret) {
        kh_value(family->members, iter) = individual;
        ret = 0;
    }
    
    // ...and founders, if applies
    if (!strcmp(individual->father_id, "0") && !strcmp(individual->mother_id, "0")) {
        ret |= family_add_founder(individual, family);
    }
    
    return ret;
}

int family_add_unknown(individual_t *individual, family_t *family) {
    assert(individual);
    assert(family);
    
    if (family_contains_individual(individual, family)) {
        // TODO return correct error code
        return -1;
    }
    
    int ret = 0;
    khiter_t iter = kh_put(family_members, family->members, strdup(individual->id), &ret);
    if (ret) {
        kh_value(family->members, iter) = individual;
        return 0;
    } else {
        return 1;
    }
}

void family_free(family_t *family) {
    assert(family);
    free(family->id);
    // Free people inside the 'members' hashtable
    for (int k = kh_begin(family->members); k < kh_end(family->members); k++) {
        if (!kh_exist(family->members, k)) { continue; }

        individual_t *ind = kh_value(family->members, k);
        individual_free(ind);
        free((void *)kh_key(family->members, k));
        kh_del(family_members, family->members, k);
    }
    
    kh_destroy(family_members, family->members);
    kh_destroy(family_members, family->founders);
    linked_list_free(family->unknown, NULL);
    
    free(family);
}

individual_t *family_contains_individual(individual_t *individual, family_t *family) {
    assert(individual->id);
    assert(family->members);
    khiter_t iter = kh_get(family_members, family->members, individual->id);
    if (iter != kh_end(family->members)) {
        return kh_value(family->members, iter);
    }
    return NULL;
}
