/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of bioinfo-libs.
 *
 * bioinfo-libs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * bioinfo-libs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bioinfo-libs. If not, see <http://www.gnu.org/licenses/>.
 */

#include "cellbase_connector.h"


char *compose_cellbase_ws_request(const char *db_host_url, const char *db_version, const char *db_species, 
                                  const char *category, const char *resource) {
    if (!db_host_url || !db_version || !db_species) {
        return NULL;
    }
    
    // URL Constants
    const char *ws_root_url = "rest/";
    const char *ws_extra_params = "?header=false&of=json";
    
    // Length of URL parts
    const int host_url_len = strlen(db_host_url);
    const int ws_root_len = strlen(ws_root_url);
    const int version_len = strlen(db_version);
    const int species_len = strlen(db_species);
    const int ws_name_len = strlen(category);
    const int resource_len = strlen(resource);
    const int ws_extra_len = strlen(ws_extra_params);
    const int result_len = host_url_len + ws_root_len + version_len + species_len + ws_name_len + resource_len + ws_extra_len + 5;
    
    char *result_url = (char*) calloc (result_len, sizeof(char));
    
    // Host URL
    strncat(result_url, db_host_url, host_url_len);
    if (result_url[host_url_len - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Root of the web service
    strncat(result_url, ws_root_url, ws_root_len);
    
    // Version
    strncat(result_url, db_version, version_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Species
    strncat(result_url, db_species, species_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Name of the web service
    strncat(result_url, category, ws_name_len);
    strncat(result_url, "/", 1);
    strncat(result_url, resource, resource_len);
    
    // Extra arguments of the web service
    strncat(result_url, ws_extra_params, ws_extra_len);
    
    return result_url;
}
