#include "array_utils.h"


typedef struct _double_index{ // Struct to store and order the values of the amplitudes preserving the index in the original array
    double value;
    size_t index;
} double_index_t;

int _compare_doubles_asc (const void *a, const void *b);

int _compare_doubles_desc (const void *a, const void *b);


int array_scalar_multiply(double *data, size_t n, double value) {
    if(data == NULL) { return -1; }
    if(n == 0) { return -2; }
    
    for(size_t i = 0; i < n; i++) {
        data[i] = data[i] * value;
    }
    return 0;
}

int array_scalar_sum(double *data, size_t n, double value) {
    if(data == NULL) { return -1; }
    if(n == 0) { return -2; }
    
    for(size_t i = 0; i < n; i++) {
        data[i] = data[i] + value;
    }
    return 0;
}


int array_sum(double *data1, const double *data2, size_t n) {
    if(data1 == NULL) { return -1; }
    if(data2 == NULL) { return -2; }
    if(n == 0) { return -3; }
    
    for(size_t i = 0; i < n; i++) {
        data1[i] = data1[i] + data2[i];
    }
    return 0;
}

int array_substract(double *data1, const double *data2, size_t n) {
    if(data1 == NULL) { return -1; }
    if(data2 == NULL) { return -2; }
    if(n == 0) { return -3; }
    
    for(size_t i = 0; i < n; i++) {
        data1[i] = data1[i] - data2[i];
    }
    return 0;
}

int array_dotproduct(double* data1, const double* data2, size_t n, double *res) {
    if(data1 == NULL) { return -1; }
    if(data2 == NULL) { return -2; }
    if(n == 0) { return -3; }
    if(res == NULL) { return -4; }
    
    *res = 0.0;
    for(size_t i = 0; i < n; i++) {
        *res += data1[i] * data2[i];
    }
    return 0;
}


int array_log(double *values, size_t n) {
    if(values == NULL) { return -1; }
    if(n == 0) { return -2; }
    
    for(size_t i = 0; i < n; i++) {
        values[i] = log(values[i]);
    }
    return 0;
}

int array_log10(double *values, size_t n) {
    if(values == NULL) { return -1; }
    if(n == 0) { return -2; }
    
    for(size_t i = 0; i < n; i++) {
        values[i] = log10(values[i]);
    }
    return 0;
}

int array_log_base(double *values, size_t n, double base) {
    if(values == NULL) { return -1; }
    if(n == 0) { return -2; }
    
    double d;
    for(size_t i = 0; i < n; i++) {
        d = log_base(values[i], base);
        if(d == -0.0) { d = 0.0; }
        values[i] = d;
    }
    return 0;
}


int array_accum(const double *values, size_t n, double *res) {
    int code = array_accum_range(values, 0, n, res);
    return (code == -4) ? -3 : code;
}

int array_accum_range(const double *values, size_t begin, size_t end, double *res) {
    if(values == NULL) { return -1; }
    if(end <= begin) { return -2; }
    if(res == NULL) { return -4; }
    
    *res = 0.0;
    for (size_t i = begin; i < end; i++) {
        if(!isnan(values[i])) {
            *res += values[i];
        }
    }
    return 0;
}

int array_order(double *values, size_t n, int asc, size_t *indices) {
    double_index_t *array_doubles_indices = (double_index_t*) malloc(n * sizeof(double_index_t));
    
    /*  Initialize struct to sort    */
    for(size_t i = 0; i < n; i++){
        array_doubles_indices[i].value = values[i];
        array_doubles_indices[i].index = i;
    }
    
    /*  Sort ASC or DESC    */
    if(asc == 1) {
        qsort(array_doubles_indices, n, sizeof(array_doubles_indices[0]), _compare_doubles_asc);
    }else {
        qsort(array_doubles_indices, n, sizeof(array_doubles_indices[0]), _compare_doubles_desc);
    }

    /*  Save sorted indices */
    for(size_t i = 0; i < n; i++){
        indices[i] = array_doubles_indices[i].index;
    }
    
    free(array_doubles_indices);
    return 0;
}

int array_ordered(const double *values, size_t n, const size_t *indices, double *ordered) {
    for(size_t i = 0; i < n; i++) {
        ordered[i] = values[indices[i]];
    }
    return 1;
}

void array_shuffle(double *values, size_t n) {
    if (n > 1) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        int usec = tv.tv_usec;
        srand48(usec);
        
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48() * (i+1));
            double t = values[j];
            values[j] = values[i];
            values[i] = t;
        }
    }
}

void array_shuffle_int(int *values, size_t n) {
    if (n > 1) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        int usec = tv.tv_usec;
        srand48(usec);
        
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48() * (i+1));
            int t = values[j];
            values[j] = values[i];
            values[i] = t;
        }
    }
}



size_t array_printf(const double *values, size_t n, char *format) {
    return array_fprintf(values, n, format, stdout);
}

size_t array_fprintf(const double *values, size_t n, char *format, FILE *file) {
    if(values == NULL) { return -1; }
    if(file == NULL) { return -4; }
    
    if(format == NULL) {
        format = "%f\n";
    }
    size_t num_chars = 0;
    size_t i;
    for(i=0; i < n; i++) {
        num_chars += fprintf(file, format, values[i]);
    }
    return num_chars;
}

void array_fread(FILE *file, double *values, size_t n) {
    size_t cont = 0;
    char line[40];
    while(cont < n && fgets(line, 40, file) != NULL) {
        values[cont] = atof(line);
        cont++;
    }
}



int _compare_doubles_asc (const void *a, const void *b) {
    double_index_t *struct_a = (double_index_t*) a;
    double_index_t *struct_b = (double_index_t*) b;

    if (struct_a->value > struct_b->value) return 1;
    else if (struct_a->value == struct_b->value) return 0;
    else return -1;
}

int _compare_doubles_desc (const void *a, const void *b) {
    double_index_t *struct_a = (double_index_t*) a;
    double_index_t *struct_b = (double_index_t*) b;

    if (struct_a->value < struct_b->value) return 1;
    else if (struct_a->value == struct_b->value) return 0;
    else return -1;
}