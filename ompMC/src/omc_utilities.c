/******************************************************************************
 ompMC - An OpenMP parallel implementation for Monte Carlo particle transport
 simulations
 
 Copyright (C) 2020 Edgardo Doerner (edoerner@fis.puc.cl)


 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
*****************************************************************************/

#include "omc_utilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Redefine printf() function due to conflicts with mex and OpenMP */
#ifdef _OPENMP
    #include <omp.h>

    #undef printf
    #define printf(...) fprintf(stderr,__VA_ARGS__)
#endif

/******************************************************************************/
/* Timing utilities. If OpenMP is enabled it calculates the wall time through 
 omp_get_wtime() function. Otherwise, it calculates CPU time through the clock() 
 function, available in time.h library. */

double omc_get_time() {

    double time_s;    // time in seconds.

#ifdef _OPENMP
    time_s = omp_get_wtime();
#else
    time_s = (double)clock()/CLOCKS_PER_SEC;
#endif

    return time_s;
}
/******************************************************************************/

/******************************************************************************/
/* A simple C/C++ class to parse input files and return requested
 key value -- https://github.com/bmaynard/iniReader */

#include <string.h>
#include <ctype.h>

/* Parse a configuration file */
void parseInputFile(char *input_file) {
    
    char buf[BUFFER_SIZE];      // support lines up to 120 characters
    
    /* Make space for the new string */
    const char *extension = INPUT_EXT;
    char *file_name = malloc(strlen(input_file) + strlen(extension) + 1);
    strcpy(file_name, input_file);
    strcat(file_name, extension); /* add the extension */
    
    FILE *fp;
    if ((fp = fopen(file_name, "r")) == NULL) {
        printf("Unable to open file: %s\n", file_name);
        exit(EXIT_FAILURE);
    }
    
    while (fgets(buf, BUFFER_SIZE , fp) != NULL) {
        /* Jumps lines labeled with #, together with only white
         spaced or empty ones. */
        if (strstr(buf, "#") || lineBlack(buf)) {
            continue;
        }
        
        strcpy(input_items[input_idx].key, strtok(buf, "=\r\n"));
        strcpy(input_items[input_idx].value, strtok(NULL, "\r\n"));
        input_idx++;
    }
    
    input_idx--;
    fclose(fp);
    
    if(verbose_flag) {
        for (int i = 0; i<input_idx; i++) {
            printf("key = %s, value = %s\n", input_items[i].key,
                   input_items[i].value);
        }
    }

    /* Cleaning */
    free(file_name);
    
    return;
}

/* Copy the value of the selected input item to the char pointer */
int getInputValue(char *dest, char *key) {
    
    /* Check to see if anything got parsed */
    if (input_idx == 0) {
        return 0;
    }
    
    for (int i = 0; i <= input_idx; i++) {
        if (strstr(input_items[i].key, key)) {
            strcpy(dest, input_items[i].value);
            return 1;
        }
    }
    
    return 0;
}

/* Returns nonzero if line is a string containing only whitespace or is empty */
int lineBlack(char *line) {
    char * ch;
    int is_blank = 1;
    
    /* Iterate through each character. */
    for (ch = line; *ch != '\0'; ++ch) {
        if (!isspace(*ch)) {
            /* Found a non-whitespace character. */
            is_blank = 0;
            break;
        }
    }
    
    return is_blank;
}

/* Remove white spaces from string str_untrimmed and saves the results in
 str_trimmed. Useful for string input values, such as file names */
 void removeSpaces(char* str_trimmed,
                  const char* str_untrimmed) {
    
    while (*str_untrimmed != '\0') {
        if(!isspace(*str_untrimmed)) {
            *str_trimmed = *str_untrimmed;
            str_trimmed++;
        }
        str_untrimmed++;
    }
    
    *str_trimmed = '\0';
    return;
}

struct inputItems input_items[INPUT_PAIRS];     // key,value pairs
int input_idx = 0;                              // number of key,value pair

/******************************************************************************/

/* Flag set by '--verbose' argument */
int verbose_flag;