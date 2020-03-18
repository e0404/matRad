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

/*******************************************************************************
* Implementation, based on the EGSnrc one, of the RANMAR random number 
* generator (RNG), proposed by Marsaglia and Zaman. 
* 
* Following the EGSnrc implementation, it uses integers to store the state of 
* the RNG and to generate the next number in the sequence. Only at the end the 
* random numbers are converted to reals, due to performance reasons. 
* 
* Before using the RNG, it is needed to initialize the RNG by a call to 
* initRandom(). 
*******************************************************************************/

#include "omc_random.h"

#include <stdlib.h>

/* Common functions and definitions */
#if defined(_MSC_VER)
	/* use __declspec(thread) instead of threadprivate to avoid 
	error C3053. More information in:
	https://stackoverflow.com/questions/12560243/using-threadprivate-directive-in-visual-studio */
	__declspec(thread) struct Random rng;
#else
	#pragma omp threadprivate(rng)
	struct Random rng;
#endif

/* Initialization function for the RANMAR random number generator (RNG) 
proposed by Marsaglia and Zaman and adapted from the EGSnrc version to be 
used in ompMC. */
void initRandom() {

    int ixx, jxx;
    
    /* Get initial seeds from input */
    char buffer[BUFF_SIZE];
    if (getInputValue(buffer, "rng seeds") != 1) {
        printf("Can not find 'rng seeds' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    sscanf(buffer, "%d %d", &ixx, &jxx);

    /* Modify jxx seed depending on OpenMP thread id */
    #ifdef _OPENMP
        jxx = jxx + omp_get_thread_num();
    #endif
    
    if (ixx <= 0 || ixx > 31328) {
        printf("Warning!, setting Marsaglia default for ixx\n");
        ixx = 1802; /* sets Marsaglia default */
    }
    if (jxx <= 0 || jxx > 31328) {
        printf("Warning!, setting Marsaglia default for jxx\n");
        jxx = 9373; /* sets Marsaglia default */
    }

    /* Save seeds to rng state struct and print information to console */
    rng.ixx = ixx;
    rng.jxx = jxx;

    printf("RNG seeds : ixx = %d, jxx = %d\n", rng.ixx, rng.jxx);
    
    int i = (rng.ixx/177 % 177) + 2;
    int j = (rng.ixx % 177) + 2;
    int k = (rng.jxx/169 % 178) + 1;
    int l = (rng.jxx % 169);
    
    int s, t, m;
    rng.urndm = malloc(97*sizeof(int));
    for (int ii = 0; ii<97; ii++) {
        s = 0;
        t = 8388608;    /* t is 2^23, half of the maximum allowed. Note that
                         only 24 bits are used */
        for (int jj=0; jj<24; jj++) {
            m = ((i*j % 179)*k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53*l + 1) % 169;
            
            if (l*m % 64 >= 32) {
                s += t;
            }
            t /=2;
        }
        rng.urndm[ii] = s;
    }
    
    rng.crndm = 362436;
    rng.cdrndm = 7654321;
    rng.cmrndm = 16777213;
    
    rng.twom24 = 1.0/16777216.0;
    
    rng.ixx = 97;
    rng.jxx = 33;
    
    /* Allocate memory for random array and set seed to start calculation of
     random numbers */
    rng.rng_array = malloc(NRANDOM*sizeof(int));
    rng.rng_seed = NRANDOM;
    
    return;
}

/* Generation function for the RANMAR random number generator (RNG) proposed 
by Marsaglia and Zaman. It generates NRANDOM floating point numbers in 
each call */
void getRandom() {
    
    int iopt;
    
    for (int i=0; i<NRANDOM; i++) {
        iopt = rng.urndm[rng.ixx - 1] - rng.urndm[rng.jxx - 1]; /* C index */
        if (iopt < 0) {
            iopt += 16777216;
        }
        
        rng.urndm[rng.ixx - 1] = iopt;

        rng.ixx -= 1;
        rng.jxx -= 1;
        if (rng.ixx == 0) {
            rng.ixx = 97;
        }
        else if (rng.jxx == 0) {
            rng.jxx = 97;
        }
        
        rng.crndm -= rng.cdrndm;
        if (rng.crndm < 0) {
            rng.crndm += rng.cmrndm;
        }
        
        iopt -= rng.crndm;
        if (iopt < 0) {
            iopt += 16777216;
        }
        rng.rng_array[i] = iopt;
    }
    
    rng.rng_seed = 0; /* index in C starts at 0 */
    
    return;
}

/* Get a single floating random number in [0,1) using the RANMAR RNG */
double setRandom() {
    
    double rnno = 0.0;
    
    if (rng.rng_seed >= NRANDOM) {
        getRandom();
    }
    
    rnno = rng.rng_array[rng.rng_seed]*rng.twom24;
    rng.rng_seed += 1;
    
    return rnno;
}

void cleanRandom() {
    
    free(rng.urndm);
    free(rng.rng_array);
    
    return;
}

/******************************************************************************/