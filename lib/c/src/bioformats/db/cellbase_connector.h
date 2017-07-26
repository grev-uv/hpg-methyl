/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CELLBASE_CONNECTOR_H
#define	CELLBASE_CONNECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Given a list of arguments, composes a URL to invoke a web service.
 * @detail Given a list of arguments, compounds a URL to invoke a web service.
 * 
 * @param db_host_url URL of the host where the database runs
 * @param db_version Version of the database to query
 * @param db_species Species whose information will be queried
 * @param category Category and subcategory of the URL (such as "genomic/variant")
 * @param resource Resource to query (such as "consequence_type")
 * 
 * @return URL composed from the given arguments, or NULL if any of the mandatory arguments is NULL
 */
char *compose_cellbase_ws_request(const char *db_host_url, const char *db_version, const char *db_species, 
                                  const char *category, const char *resource);

#ifdef __cplusplus
}
#endif

#endif
