#ifndef OMC_UTILITIES_H
#define OMC_UTILITIES_H
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

/******************************************************************************/
/* Timing utilities. If OpenMP is enabled it calculates the wall time through 
 omp_get_wtime() function. Otherwise, it calculates CPU time through the clock() 
 function, available in time.h library. */

double omc_get_time();
/******************************************************************************/

/******************************************************************************/
/* A simple C/C++ class to parse input files and return requested key value 
https://github.com/bmaynard/iniReader */

#define BUFFER_SIZE 256
#define INPUT_PAIRS 80
#define INPUT_EXT ".inp"  // extension of input files

/* Parse a configuration file */
extern void parseInputFile(char *file_name);

/* Copy the value of the selected input item to the char pointer */
extern int getInputValue(char *dest, char *key);

/* Returns nonzero if line is a string containing only whitespace or is empty */
extern int lineBlack(char *line);

/* Remove white spaces from string str_untrimmed and saves the results in
 str_trimmed. Useful for string input values, such as file names */
extern void removeSpaces(char* str_trimmed, const char* str_untrimmed);

struct inputItems {
    char key[BUFFER_SIZE];
    char value[BUFFER_SIZE];
};

extern struct inputItems input_items[];     // key,value pairs
extern int input_idx;                       // number of key,value pair

/******************************************************************************/

/* Flag set by '--verbose' argument */
extern int verbose_flag;

#endif