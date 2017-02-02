#ifndef HTTP_UTILS_H
#define HTTP_UTILS_H

#include <stdlib.h>
#include <string.h>

#include <curl/curl.h>

/**
 * Initialize the environment, setting whether SSL should be active or not.
 * 
 * @param ssl If the authentication will be performed via SSL
 */
int init_http_environment(int ssl);

/**
 * Request a URL with the provided parameters using HTTP GET method. A callback function
 * can be provided for its invocation when the response is served.
 * 
 * @param url URL to request
 * @param params GET parameters of the request
 * @param write_function Function for managing the response contents
 */
int http_get(char *url, char **params, char **params_values, int num_params, size_t (*write_function) (char*, size_t, size_t, void*), void *buffer);

/**
 * Request a URL with the provided parameters using HTTP POST method. A callback function
 * can be provided for its invocation when the response is served.
 * 
 * @param url URL to request
 * @param params POST parameters of the request
 * @param write_function Function for managing the response contents
 */
int http_post(char *url, char **params, char **params_values, int num_params, size_t (*write_function) (char*, size_t, size_t, void*));

/**
 * Should an error occur, return the string message associated to a numerical code.
 * 
 * @param err_code Numerical error code
 * 
 * @return Message related to the provided error code
 */
const char *get_last_http_error(int err_code);

#endif
