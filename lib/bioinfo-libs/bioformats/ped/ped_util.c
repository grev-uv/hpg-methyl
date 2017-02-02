#include "ped_util.h"

enum Condition get_condition_from_phenotype(float phenotype) {
    if (phenotype - 1.0 < 1e-6) {
        return UNAFFECTED;
    } else if (phenotype - 2.0 < 1e-6) {
        return AFFECTED;
    } else {
        return UNKNOWN_CONDITION;
    }
}
