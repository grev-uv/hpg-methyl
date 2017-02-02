#include "file_utils.h"

/* **************************************************************
 *    		Functions implementations  			*
 * **************************************************************/
 
void *mmap_file(size_t *len, const char *filename) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        LOG_FATAL_F("Error opening file: %s\n", filename);
    }
    
    struct stat st[1];
    if (fstat(fd, st)) {
        LOG_FATAL_F("Error while getting file information: %s\n", filename);
    }
    *len = (size_t) st->st_size;

    if (!*len) {
        close(fd);
        return NULL;
    }

    void *map = mmap(NULL, *len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (MAP_FAILED == map) {
        LOG_FATAL_F("mmap failed for %s\n", filename);
    }
    close(fd);
    
    return map;
}


char *fgets_no_ln(char *s, int n, FILE *f) {
    int c = 0;
    char *cs;

    cs = s;
    while (--n > 0 && (c = getc(f)) != EOF) {
        if (c == '\n') {
            break;
        }
        *cs++ = c;
    }
    if (c == EOF && cs == s)
        return NULL;
    *cs++ = '\0';
    return s;
}

int exists(const char * filename) {
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return 1;
    }
    return 0;
}

unsigned long count_lines(const char *filename) {
    unsigned long count = 0;
    char *buffer = (char*) malloc(1024 * 4);
    FILE *file = fopen(filename, "r");
    while (fgets(buffer, 1024 * 4, file) != NULL) {
        count++;
    }
    fclose(file);
    return count;
}

int create_directory(const char* path) {
    return mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

int delete_directory(const char* path) {
    return rmdir(path);
}

int delete_files_by_extension(const char* dir_path, const char* extension) {
    DIR *dp;
    struct dirent *ep;
    char *file_path;

    dp = opendir(dir_path);

    if (dp == NULL) {
        return 1;
    }

    while ((ep = readdir(dp)) != NULL) {
        if (ends_with(ep->d_name, extension)) {
            file_path = (char*) calloc(strlen(dir_path) + strlen(ep->d_name) + 2, sizeof(char));
            sprintf(file_path, "%s/%s", dir_path, ep->d_name);
//             strncat(file_path, dir_path, strlen(dir_path));
//             strncat(file_path, ep->d_name, strlen(ep->d_name));
            remove(file_path);
            free(file_path);
        }
    }

    return closedir(dp);
}

char** parse_conf_file(char *filename) {
    unsigned long num_lines = count_lines(filename);
    char ** option_line = (char **)malloc(2 * num_lines * sizeof(char *));
    const char delimiters[] = "=+";

    FILE *fd = fopen(filename, "r");

    int option_count = 0;
    char line[MAX_LENGTH_CONFIG_LINE];
    while (fgets(line, MAX_LENGTH_CONFIG_LINE, fd) != NULL) {

        if (line[0] != '#' && line[0] != '\n') {
            chomp(line);
            strcpy(line, str_replace(line, (const char*)" ", (const char*)"", MAX_LENGTH_CONFIG_LINE));

            option_line[option_count] = (char *)malloc((strlen(line) + 2) * sizeof(char) + 1);
            strcpy(option_line[option_count], "--");
            strcat(option_line[option_count], line);

            char *option_line_copy;

            option_line_copy = strdup(option_line[option_count]);
            sscanf(strtok(option_line_copy, delimiters), "%s", option_line[option_count]);
            option_count++;

            option_line[option_count] = (char *)malloc((strlen(line) + 2) * sizeof(char) + 1);
            sscanf(strtok(NULL, delimiters), "%s", option_line[option_count]);
            option_count++;
        }
    }
    fclose(fd);

    return option_line;
}

int parse_conf_file2(char **argvs, char *filename) {
    const char delimiters[] = "=";
    int linenum = 0;

    FILE *fd = fopen(filename, "r");
    char line[MAX_LENGTH_CONFIG_LINE];
    char *tokens;
    while (fgets(line, MAX_LENGTH_CONFIG_LINE, fd) != NULL) {
        if (line[0] == '#')
            continue;
        chomp(line);
        char *line_copy = strdup(line);
        tokens = strtok(line_copy, delimiters);
        argvs[linenum] = (char*)malloc(strlen((const char*)tokens) * sizeof(char));
        strcpy(argvs[linenum], tokens);
        tokens = strtok(NULL, delimiters);
        argvs[linenum+1] = (char*)malloc(strlen((const char*)tokens) * sizeof(char));
        strcpy(argvs[linenum+1], tokens);
        linenum += 2;
    }
    fclose(fd);

    return linenum;
}

char* get_filename_from_path(char* path, char* filename_p) {
    int i = 0;
    char* char_p = strrchr(path, '/');

    if (char_p == NULL) {  //the path has no subdirectories
        strcpy(filename_p, path);
    } else {
        while (*char_p != '\0') {
            char_p++;
            filename_p[i++] = *char_p;
        }
    }

    return filename_p;
}
