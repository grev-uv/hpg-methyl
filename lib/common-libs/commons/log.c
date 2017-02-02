#include "log.h"

void init_log() {
    log_level = LOG_INFO_LEVEL;
    log_verbose = 0;
    log_file = fopen("output.log", "w");
}

void init_log_custom(int level, int verbose, char* filename, char *mode) {
    log_level = level;
    log_verbose = verbose;
    log_file = fopen(filename, mode);
}

void stop_log() {
    fclose(log_file);
}

void print_log_message(int level, char *log_level_word, char *filename, int num_line, const char *func, char *msg) {
    if((level >= log_level) && (log_verbose || log_file)) {
        time_t rawtime;
        time(&rawtime);
        char* str_time = ctime(&rawtime);
        chomp(str_time);

        // if 'verbose' logs are printed in stderr
        if(log_verbose) {
            fprintf(stderr, "%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word, filename, num_line, func, msg);
        }

        // if 'log_file' has been set up then logs are printed
        // logs are ALWAYS printed in log_file independently of 'verbose' parameter
        if(log_file) {
            fprintf(log_file, "%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word, filename, num_line, func, msg);
        }

    }
}

void print_log_message_with_format(int level, char *log_level_word, char *filename, int num_line, const char *func, char *msg_format, ...) {
    if ((level >= log_level) && (log_verbose || log_file != NULL)) {
        time_t rawtime;
        time(&rawtime);
        char* str_time = ctime(&rawtime);
        chomp(str_time);

        va_list args;
        va_start(args, msg_format);
        
        // if 'verbose' logs are printed in stdout
        if (log_verbose) {
            fprintf(stderr, "%s\t%s\t%s [%i] in %s(): ", str_time, log_level_word, filename, num_line, func);

            if (args != NULL) {
            	vfprintf(stderr, msg_format, args);
            } else {
            	fprintf(stderr, "No arguments provided\n");
            }
        }

        // if 'log_file' has been set up then logs are printed
        // logs are ALWAYS printed in log_file independently of 'verbose'       
        if (log_file) {
            va_list argsf;
            memcpy(argsf, args, sizeof(args));
            va_start(argsf, msg_format);
            
            fprintf(log_file, "%s\t%s\t%s [%i] in %s(): ", str_time, log_level_word, filename, num_line, func);

            if (args != NULL) {
                vfprintf(log_file, msg_format, argsf);
            } else {
                fprintf(stderr, "No arguments provided\n");
            }

            va_end(argsf);
        }
        va_end(args);
    }
}
