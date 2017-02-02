#ifndef BAM_STATS_REPORT_H
#define BAM_STATS_REPORT_H

/*
 * bam_stats_report.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

//#include "argtable2.h"
//#include "libconfig.h"
//#include "commons/log.h"
//#include "commons/system_utils.h"
//#include "commons/file_utils.h"

//#include "commons/workflow_scheduler.h"
//#include "containers/array_list.h"
//#include "bioformats/bam-sam/bam_file.h"
//#include "bioformats/features/region/region_table.h"

#include "bam_stats.h"

//------------------------------------------------------------------------

#define MAX_SERIES_GRAPH	10
#define MAX_LINE_GRAPH  	1024

//------------------------------------------------------------------------


/**
* @brief Structure for report graphs parameters
* 
* Structure containing parameters for graphics representation of report graphs
*/
typedef struct report_graph {
    int x_autoscale;					/**< Autoscale flag for X axis. */
    int x_start;					/**< X axis start coordinate. */
    int x_end;						/**< X axis end coordinate. */
    int y_autoscale;					/**< Autoscale flag for Y axis. */
    int y_start;					/**< Y axis start coordinate. */
    int y_end;						/**< Y axis end coordinate. */
    int lmargin;					/**< Left margin. */
    int rmargin;					/**< Right margin. */	
    int tmargin;					/**< Top margin. */
    int bmargin;					/**< Bottom margin. */
    int num_y_columns;					/**< Number of columns in the Y axis. */
    int x_column;					/**< Number of column in data file with X axis data. */
    int y_columns[MAX_SERIES_GRAPH];			/**< Numbers of columns in data file that are represented in Y axis. */
    char y_titles[MAX_SERIES_GRAPH][MAX_LINE_GRAPH];			/**< Titles of series represented in Y axis. */
    char title[MAX_LINE_GRAPH];					/**< Title of the graphic. */
    char xlabel[MAX_LINE_GRAPH];					/**< Label of X axis. */	
    char ylabel[MAX_LINE_GRAPH];					/**< Label of Y axis. */
    char type[MAX_LINE_GRAPH]; 					/**< Type of graph. */
} report_graph_t;


//------------------------------------------------------------------------

void report_bam_stats_output(char *bam_filename, char *outdir, 
			     void *db, bam_stats_output_t *stats);

//------------------------------------------------------------------------

#endif // BAM_STATS_REPORT_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
