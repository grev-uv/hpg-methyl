#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "log.h"
#include "string_utils.h"

#define MAX_LENGTH_CONFIG_LINE 		256
#define MAX_FILENAME_LENGTH		256
#define MAX_FULL_PATH_LENGTH		2048

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Launch option
* 
* Structure containing an option for launching a program
*/
typedef struct launch_option {
    char option_name[MAX_LENGTH_CONFIG_LINE + 2];	/**< Name of the option. */
    char option_value[MAX_LENGTH_CONFIG_LINE];		/**< Value of the option. */
} launch_option_t;

/* **************************************
 *  		Functions		*
 * *************************************/

void *mmap_file(size_t *len, const char *filename);
 
/**
*  @brief Returns n characters from a file with no break line characters
*  @param s pointer to the char array to store the read characters
*  @param n number of characters to read
*  @param f file descriptor 
*  @return pointer to the char array with the stored characters
*  
*  Returns n characters from a file with no break line characters
*/
char* fgets_no_ln(char *s, int n, FILE *f);

/**
*  @brief Determines if a given path exists
*  @param path pointer to the path to locate
*  @return 1: the path does exist, 0: the path does not exist
*  
*  Determines if a given path exists (1: exists, 0: does not exist)
*/
int exists(const char *path);

/**
*  @brief Determines if a given path is a file
*  @param path pointer to the path to locate
*  @return 1: the path is a file, 0: the path is not a file
*  
*  Determines if a given path is a file (1: file, 0: not file)
*/
int is_file(const char *path);

/**
*  @brief Determines if a given path is a directory
*  @param path pointer to the path to locate
*  @return 1: the path is a directory, 0: the path is not a directory
*  
*  Determines if a given path is a directory (1: directory, 0: not directory)
*/
int is_directory(const char *path);

/**
*  @brief Counts lines from a file
*  @param filename pointer to the filename for line count processing
*  @return number of lines counted
*  
*  Count lines from a given file 
*/
unsigned long count_lines(const char *filename);

/**
*  @brief Copies the content from a char array to another char array
*  @param dest pointer to the destination char array
*  @param ori pointer to the origin char array
*  @return sucess value (1: copy ok, 0: copy not ok)
*  
*  Copy the content from a char array to another char array, 
*  the origin content is preserved
*/
int copy(const char *dest, const char *src);

/**
*  @brief Moves the content from a char array to another char array
*  @param dest pointer to the destination char array
*  @param ori pointer to the origin char array
*  @return sucess value (1: copy ok, 0: copy not ok)
*  
*  Move the content from a char array to another char array, 
*  the origin content is removed
*/
int move(const char *dest, const char *src);

/**
*  @brief Changes file timestamp (touch command)
*  @param path pointer to the path to touch
*  @return sucess value (1: ok, 0: not ok)
*  
*  Changes file timestamp (touch command)
*/
int touch(const char *path);

/**
*  @brief Creates a directory
*  @param path pointer to the directory path to create
*  @return sucess value (1: created, 0: not created)
*  
*  Creates a directory with a given path
*/
int create_directory(const char *path);

/**
*  @brief Deletes a directory
*  @param path pointer to the directory path to delete
*  @return sucess value (1: deleted, 0: not deleted)
*  
*  Deletes a directory with a given path
*/
int delete_directory(const char *path);

/**
*  @brief Deletes files from a directory with a given extension
*  @param path pointer to the directory from which files will be deleted
*  @param extension pointer to the file extension to delete
*  @return sucess value (1: deleted, 0: not deleted)
*  
*  Deletes files from a directory with a given extension
*/
int delete_files_by_extension(const char *dir_path, const char *extension);

/**
*  @brief Obtains stored size in a given path
*  @param path pointer to the directory to obtain the size
*  @return stored size
*  
*  Obtains stored size in a given path
*/
unsigned long size(const char *path);

/**
*  @brief Parses options from a configuration file
*  @param filename pointer to the configuration file
*  @return pointers to the parsed options
*  
*  Parses options from a configuration file that are written in option-value per line format
*/
char** parse_conf_file(char *filename);

/**
*  @brief Parses options from a configuration file
*  @param filename pointer to the configuration file
*  @return pointers to the parsed options
*  
*  Parses options from a configuration file that are written in option-value per line format
*/
int parse_conf_file2(char **argvs, char *filename);

/**
*  @brief Isolates the file name from a full path name
*  @param path pointer to the full path name
*  @param filename_p pointer to the file name
*  @return pointer to the file name
*  
*  Given a path with a filename at the end, gets only the filename
*/
char* get_filename_from_path(char* path, char* filename_p);

#endif	/*  FILE_UTILS_H   */
