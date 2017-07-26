
#include <stdio.h>

#include "../log.h"

int main(void) {
	printf("Testing  combinations of log_level and verbose\n");
	printf("==============================================\n");
	LOG_LEVEL(1);
	LOG_VERBOSE(1);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");

	LOG_LEVEL(3);
	LOG_VERBOSE(1);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");

	LOG_LEVEL(1);
	LOG_VERBOSE(0);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");
	printf("\n\n");


	printf("Testing  LOG_IF and verbose\n");
	printf("==============================================\n");
	LOG_LEVEL(LOG_DEBUG_LEVEL);
	LOG_VERBOSE(1);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_IF(LOG_DEBUG_LEVEL, (1>0), "debug message");
	LOG_IF(LOG_INFO_LEVEL, (1>0), "info message");
	LOG_IF(LOG_WARN_LEVEL, (1>0), "warn message");
	LOG_IF(LOG_ERROR_LEVEL, (1>0), "error message"
	);
	printf("\n\n");

	printf("Testing  combinations of log_level and verbose and writing to\
/tmp/test_log.log\n");
	printf("==============================================\n");
	LOG_LEVEL(1);
	LOG_VERBOSE(1);
	LOG_FILE("/tmp/test_log.log", "w");
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");

	LOG_LEVEL(3);
	LOG_VERBOSE(1);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");

	LOG_LEVEL(1);
	LOG_VERBOSE(0);
	printf("Testing logs with log_level=%i and verbose=%i\n", log_level,
log_verbose);
	LOG_DEBUG("debug message");
	LOG_INFO("info message");
	LOG_WARN("warn message");
	LOG_ERROR("error message");
	printf("\n\n");

	return 1;
}