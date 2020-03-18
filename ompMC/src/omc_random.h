#ifndef OMC_RANDOM_H
#define OMC_RANDOM_H
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
* random numbers are converted to floating point numbers, due to performance 
* reasons. 
* 
* Before using the RNG, it is needed to initialize the RNG by a call to 
* initRandom(). 
*******************************************************************************/

#define NRANDOM 128     // number of random numbers generated in each call 
                        // to setRandom().
#define BUFF_SIZE 256

struct Random {
    int crndm;
    int cdrndm;
    int cmrndm;
    int ixx;
    int jxx;
    int rng_seed;
    int *urndm;
    int *rng_array;
    double twom24;
};

#if defined(_MSC_VER)
	/* use __declspec(thread) instead of threadprivate to avoid
	error C3053. More information in:
	https://stackoverflow.com/questions/12560243/using-threadprivate-directive-in-visual-studio */
	extern __declspec(thread) struct Random rng;
#else
	extern struct Random rng;
	#pragma omp threadprivate(rng)
#endif

/* Initialization function for the RANMAR random number generator (RNG) 
proposed by Marsaglia and Zaman and adapted from the EGSnrc version to be 
used in ompMC. */
extern void initRandom(void);

/* Generation function for the RANMAR random number generator (RNG) proposed 
by Marsaglia and Zaman. It generates NRANDOM floating point numbers in 
each call */
extern void getRandom(void);

/* Get a single floating random number in [0,1) using the RANMAR RNG */
extern double setRandom(void);

extern void cleanRandom(void);

/******************************************************************************/

#endif