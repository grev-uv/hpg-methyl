#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../file_utils.h"

int main(void) {
	printf("Testing  combinations of log_level and verbose\n");
	printf("==============================================\n");
	char *filename = (char*)malloc(70*sizeof(char));
	strcpy(filename, "/home/imedina/appl/bioinfo-c/commons/file_utils.h");
	unsigned long num_lines = count_lines(filename);
	printf("%lu", num_lines);

	return 1;
}