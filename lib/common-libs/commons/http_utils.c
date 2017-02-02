#include "http_utils.h"


int init_http_environment(int ssl) {
    if (ssl) {
        return curl_global_init(CURL_GLOBAL_SSL);
    } else {
        return curl_global_init(CURL_GLOBAL_NOTHING);
    }
}

int http_get(char *url, char **params, char **params_values, int num_params, size_t (*write_function) (char*, size_t, size_t, void*), void *buffer) {
    CURL *curl;
    CURLcode ret_code = CURLE_OK;

    struct curl_httppost *formpost = NULL;
    char *aux_url;
    
    if (params != NULL && params_values != NULL) {
        // Add request parameters to the URL
        for (int i = 0; i < num_params; i++) {
            int param_len;
            if (i > 0) { param_len = strlen(params[i]) + strlen(params_values[i]) + 1; }
            else       { param_len = strlen(params[i]) + strlen(params_values[i]) + 2; }
            aux_url = (char*) realloc (url, (strlen(url) + param_len + 1) * sizeof(char));
            if (aux_url) {
                url = aux_url;
                if (i == 0) { strncat(url, "?", 1); }
                strncat(url, params[i], strlen(params[i]));
                strncat(url, "=", 1);
                strncat(url, params_values[i], strlen(params_values[i]));
            }
        }
    }
    
    curl = curl_easy_init();
    
    if(curl) {
//         if (debug) { curl_easy_setopt(curl, CURLOPT_VERBOSE, 1); }
        
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_function);
	// TODO CURLOPT_WRITEDATA
	
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, buffer);

	
        // Option for avoiding a segfault caused by a signal in libcurl
        // http://stackoverflow.com/questions/9191668/error-longjmp-causes-uninitialized-stack-frame
        curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1);
        
        ret_code = curl_easy_perform(curl);
        if (ret_code != 0) { fprintf(stderr, "ret_code = %s\n", curl_easy_strerror(ret_code)); }

        curl_easy_cleanup(curl);
        curl_formfree(formpost);
    } else {
        ret_code = CURLE_FAILED_INIT;
    }
    
    return ret_code;
}


int http_post(char *url, char **params, char **params_values, int num_params, size_t (*write_function) (char*, size_t, size_t, void*)) {
    CURL *curl;
    CURLcode ret_code = CURLE_OK;

    struct curl_httppost *formpost = NULL;
    struct curl_httppost *lastptr = NULL;
    
    // Set request parameters
    for (int i = 0; i < num_params; i++) {
        curl_formadd(&formpost,
                    &lastptr,
                    CURLFORM_COPYNAME, params[i],
                    CURLFORM_PTRCONTENTS, params_values[i],
                    CURLFORM_END);
    }
    
    curl = curl_easy_init();
    
    if(curl) {
//         if (debug) { curl_easy_setopt(curl, CURLOPT_VERBOSE, 1); }
        
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_HTTPPOST, formpost);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_function);
        // Option for avoiding a segfault caused by a signal in libcurl
        // http://stackoverflow.com/questions/9191668/error-longjmp-causes-uninitialized-stack-frame
        curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1);
        
        ret_code = curl_easy_perform(curl);
        if (ret_code != 0) { fprintf(stderr, "ret_code = %s\n", curl_easy_strerror(ret_code)); }

        curl_easy_cleanup(curl);
        curl_formfree(formpost);
    } else {
        ret_code = CURLE_FAILED_INIT;
    }
    
    return ret_code;
}

const char *get_last_http_error(int err_code) {
    return curl_easy_strerror(err_code);
}
