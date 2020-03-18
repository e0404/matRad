/******************************************************************************
 ompMC - An OpenMP parallel implementation for Monte Carlo particle transport
 simulations
 
 Copyright (C) 2018 Edgardo Doerner (edoerner@fis.puc.cl)


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

/******************************************************************************
 omc_matrad - An ompMC user code to calculate deposited dose on voxelized 
 geometries to be used with the matRad treatment planning system.  
*****************************************************************************/

/******************************************************************************
/* Definitions needed if source file compiled with mex. This macro must be 
/* enabled during compilation time.
*****************************************************************************/
#include <mex.h>

/* Redefine printf() function due to conflicts with mex and OpenMP */
#ifdef _OPENMP
    #include <omp.h>

    #undef printf
    #define printf(...) fprintf(stderr,__VA_ARGS__)
#endif

#define exit(EXIT_FAILURE) mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:invalid","Abort.");

#include "omc_utilities.h"
#include "ompmc.h"
#include "omc_random.h"

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Variables needed to parse inputs from matRad */
const mxArray *cubeRho;
const mxArray *cubeMatIx;
const mxArray *mcGeo;
const mxArray *mcSrc;
const mxArray *mcOpt;

/* Function used to parse input from matRad */
void parseInput(int nrhs, const mxArray *prhs[]) {

    mxArray *tmp_fieldpointer;
    char *tmp;

    cubeRho = prhs[0];
    cubeMatIx = prhs[1];
    mcGeo = prhs[2];
    mcSrc = prhs[3];
    mcOpt = prhs[4];

    /* Check data type of input arguments */
    if (!(mxIsDouble(cubeRho))){
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:inputNotDouble",
                "Input argument must be of type double.");
    }    
    if (mxGetNumberOfDimensions(cubeRho) != 3){
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:inputNot3D",
                "Input argument 1 must be a three-dimensional cube\n");
    }
    if (!mxIsInt32(cubeMatIx)) {
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:inputNotInt32","The density cube must be a 32 bit integer array!");
    }
    if(!mxIsStruct(mcGeo)) {
        mexErrMsgIdAndTxt( "MATLAB:phonebook:inputNotStruct",
                "Input 3 must be a mcGeo Structure.");
    }
    if(!mxIsStruct(mcSrc)) {
        mexErrMsgIdAndTxt( "MATLAB:phonebook:inputNotStruct",
                "Input 4 must be a mcSrc Structure.");
    }

    /* Parse Monte Carlo options and create input items structure */
    tmp_fieldpointer = mxGetField(mcOpt,0,"verbose");
    verbose_flag = mxGetLogicals(tmp_fieldpointer)[0];

    mxArray* tmp2;
    int status;
    int nInput = 0;
    
    sprintf(input_items[nInput].key,"ncase");
    tmp_fieldpointer = mxGetField(mcOpt,0,"nHistories");
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    nInput++;
    sprintf(input_items[nInput].key,"nbatch");
    tmp_fieldpointer = mxGetField(mcOpt,0,"nBatches");
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    nInput++; /* Get splitting factor */
    sprintf(input_items[nInput].key,"nsplit");
    tmp_fieldpointer = mxGetField(mcOpt,0,"nSplit");
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }

    nInput++;
    sprintf(input_items[nInput].key,"spectrum file");
    tmp_fieldpointer = mxGetField(mcOpt,0,"spectrumFile");
    tmp = mxArrayToString(tmp_fieldpointer);
    strcpy(input_items[nInput].value,tmp);
    
    
    nInput++;
    sprintf(input_items[nInput].key,"mono energy");
    tmp_fieldpointer = mxGetField(mcOpt,0,"monoEnergy");
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");  

    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    nInput++;
    sprintf(input_items[nInput].key,"charge");
    tmp_fieldpointer = mxGetField(mcOpt,0,"charge");    
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }

    nInput++;
    sprintf(input_items[nInput].key,"global ecut");
    tmp_fieldpointer = mxGetField(mcOpt,0,"global_ecut");    
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }

    nInput++;
    sprintf(input_items[nInput].key,"global pcut");
    tmp_fieldpointer = mxGetField(mcOpt,0,"global_pcut");    
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    nInput++;
    sprintf(input_items[nInput].key,"rng seeds");
    tmp_fieldpointer = mxGetField(mcOpt,0,"randomSeeds");    
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    nInput++;
    sprintf(input_items[nInput].key,"pegs file");
    tmp_fieldpointer = mxGetField(mcOpt,0,"pegsFile");    
    tmp = mxArrayToString(tmp_fieldpointer);
    strcpy(input_items[nInput].value,tmp);
    
    nInput++;
    sprintf(input_items[nInput].key,"pgs4form file");
    tmp_fieldpointer = mxGetField(mcOpt,0,"pgs4formFile");    
    tmp = mxArrayToString(tmp_fieldpointer);
    strcpy(input_items[nInput].value,tmp);
    
    nInput++;
    sprintf(input_items[nInput].key,"data folder");
    tmp_fieldpointer = mxGetField(mcOpt,0,"dataFolder");    
    tmp = mxArrayToString(tmp_fieldpointer);
    strcpy(input_items[nInput].value,tmp);
    
    nInput++;
    sprintf(input_items[nInput].key,"output folder");
    tmp_fieldpointer = mxGetField(mcOpt,0,"outputFolder");    
    tmp = mxArrayToString(tmp_fieldpointer);
    strcpy(input_items[nInput].value,tmp);

    nInput++;
    sprintf(input_items[nInput].key,"relative dose threshold");
    tmp_fieldpointer = mxGetField(mcOpt,0,"relDoseThreshold");
    status = mexCallMATLAB(1, &tmp2, 1,  &tmp_fieldpointer, "num2str");    
    if (status != 0)
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Call to num2str not successful");
    else
    {
        tmp = mxArrayToString(tmp2);        
        strcpy(input_items[nInput].value,tmp);
    }
    
    input_idx = nInput;
    
    mexPrintf("Input Options:\n");
    for (int iInput = 0; iInput < nInput; iInput++)
        mexPrintf("%s: %s\n",input_items[iInput].key,input_items[iInput].value);
    
    if (verbose_flag)
        mexPrintf("ompMC output Option: Verbose flag is set!\n");
            
    return;
}

/******************************************************************************/
/* Geometry definitions */
struct Geom {
    int *med_indices;           // index of the media in each voxel
    double *med_densities;      // density of the medium in each voxel
    
    int isize;                  // number of voxels on each direction
    int jsize;
    int ksize;
    
    double *xbounds;            // boundaries of voxels on each direction
    double *ybounds;
    double *zbounds;
};
struct Geom geometry;

void initPhantom() {
    
    /* Get phantom information from matRad */
    unsigned int nfields;
    int ngeostructfields;
    mwSize ndim, nmaterials;
    mwSize *materialdim;
    mxArray *tmp_fieldpointer;

    ngeostructfields = mxGetNumberOfFields(mcGeo);

    /* Get number of media and media names. This info is saved in media struct */
    tmp_fieldpointer = mxGetField(mcGeo,0,"material");
    materialdim = mxGetDimensions(tmp_fieldpointer);
    nmaterials = materialdim[0];    
    media.nmed = nmaterials;
    
    mwIndex tmpSubs[2];
    char *tmp;
    for (int iMat = 0; iMat < nmaterials; iMat++) {
        tmpSubs[0] = iMat;
        tmpSubs[1] = 0;
        mwSize linIx = mxCalcSingleSubscript(tmp_fieldpointer,2,tmpSubs);
        mxArray* tmpCellPointer = mxGetCell(tmp_fieldpointer,linIx);
        
        tmp = mxArrayToString(tmpCellPointer);
        if (tmp)
        {
            strcpy(media.med_names[iMat],tmp);
        }
        else
        {
            mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:Error","Material string could not be read!");
        }
    }

    /* Get boundaries, density and material index for each voxel */
    const mwSize *cubeDim = mxGetDimensions(cubeRho);        
    mwSize nCubeElements = cubeDim[0]*cubeDim[1]*cubeDim[2];
    
    geometry.isize = cubeDim[0];
    geometry.jsize = cubeDim[1];
    geometry.ksize = cubeDim[2];
    
    tmp_fieldpointer = mxGetField(mcGeo,0,"xBounds");    
    geometry.xbounds = mxGetPr(tmp_fieldpointer);
    tmp_fieldpointer = mxGetField(mcGeo,0,"yBounds");    
    geometry.ybounds = mxGetPr(tmp_fieldpointer);
    tmp_fieldpointer = mxGetField(mcGeo,0,"zBounds");    
    geometry.zbounds = mxGetPr(tmp_fieldpointer);
    
    geometry.med_densities = mxGetPr(cubeRho);
    
    geometry.med_indices = (int*)mxGetPr(cubeMatIx);

    /* Summary with geometry information */
    mexPrintf("Number of media in phantom : %d\n", media.nmed);
    mexPrintf("Media names: ");
    for (int i=0; i<media.nmed; i++) {
        mexPrintf("%s, ", media.med_names[i]);
    }
    mexPrintf("\n");
    mexPrintf("Number of voxels on each direction (X,Y,Z) : (%d, %d, %d)\n",
           geometry.isize, geometry.jsize, geometry.ksize);
    mexPrintf("Minimum and maximum boundaries on each direction : \n");
    mexPrintf("\tX (cm) : %lf, %lf\n",
           geometry.xbounds[0], geometry.xbounds[geometry.isize]);
    mexPrintf("\tY (cm) : %lf, %lf\n",
           geometry.ybounds[0], geometry.ybounds[geometry.jsize]);
    mexPrintf("\tZ (cm) : %lf, %lf\n",
           geometry.zbounds[0], geometry.zbounds[geometry.ksize]);
    
    return;
}

void cleanPhantom() {
    
    /* The memory inside geometry structure is shared with Matlab, therefore 
    it is not freed here */
    
    return;
}

void howfar(int *idisc, int *irnew, double *ustep) {
    
    int np = stack.np;
    int irl = stack.ir[np];
    double dist = 0.0;
    
    if (stack.ir[np] == 0) {
        /* The particle is outside the geometry, terminate history */
        *idisc = 1;
        return;
    }
    
    /* If here, the particle is in the geometry, do transport checks */
    int ijmax = geometry.isize*geometry.jsize;
    int imax = geometry.isize;
    
    /* First we need to decode the region number of the particle in terms of
     the region indices in each direction */
    int irx = (irl - 1)%imax;
    int irz = (irl - 1 - irx)/ijmax;
    int iry = ((irl - 1 - irx) - irz*ijmax)/imax;
    
    /* Check in z-direction */
    if (stack.w[np] > 0.0) {
        /* Going towards outer plane */
        dist = (geometry.zbounds[irz+1] - stack.z[np])/stack.w[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (irz != (geometry.ksize - 1)) {
                *irnew = irl + ijmax;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }
    
    else if (stack.w[np] < 0.0) {
        /* Going towards inner plane */
        dist = -(stack.z[np] - geometry.zbounds[irz])/stack.w[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (irz != 0) {
                *irnew = irl - ijmax;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }

    /* Check in x-direction */
    if (stack.u[np] > 0.0) {
        /* Going towards positive plane */
        dist = (geometry.xbounds[irx+1] - stack.x[np])/stack.u[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (irx != (geometry.isize - 1)) {
                *irnew = irl + 1;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }
    
    else if (stack.u[np] < 0.0) {
        /* Going towards negative plane */
        dist = -(stack.x[np] - geometry.xbounds[irx])/stack.u[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (irx != 0) {
                *irnew = irl - 1;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }
    
    /* Check in y-direction */
    if (stack.v[np] > 0.0) {
        /* Going towards positive plane */
        dist = (geometry.ybounds[iry+1] - stack.y[np])/stack.v[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (iry != (geometry.jsize - 1)) {
                *irnew = irl + imax;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }
    
    else if (stack.v[np] < 0.0) {
        /* Going towards negative plane */
        dist = -(stack.y[np] - geometry.ybounds[iry])/stack.v[np];
        if (dist < *ustep) {
            *ustep = dist;
            if (iry != 0) {
                *irnew = irl - imax;
            }
            else {
                *irnew = 0; /* leaving geometry */
            }
        }
    }
    
    return;
}

double hownear(void) {
    
    int np = stack.np;
    int irl = stack.ir[np];
    double tperp = 1.0E10;  /* perpendicular distance to closest boundary */
    
    if (irl == 0) {
        /* Particle exiting geometry */
        tperp = 0.0;
    }
    else {
        /* In the geometry, do transport checks */
        int ijmax = geometry.isize*geometry.jsize;
        int imax = geometry.isize;
        
        /* First we need to decode the region number of the particle in terms
         of the region indices in each direction */
        int irx = (irl - 1)%imax;
        int irz = (irl - 1 - irx)/ijmax;
        int iry = ((irl - 1 - irx) - irz*ijmax)/imax;
        
        /* Check in x-direction */
        tperp = fmin(tperp, geometry.xbounds[irx+1] - stack.x[np]);
        tperp = fmin(tperp, stack.x[np] - geometry.xbounds[irx]);
        
        /* Check in y-direction */
        tperp = fmin(tperp, geometry.ybounds[iry+1] - stack.y[np]);
        tperp = fmin(tperp, stack.y[np] - geometry.ybounds[iry]);
        
        /* Check in z-direction */
        tperp = fmin(tperp, geometry.zbounds[irz+1] - stack.z[np]);
        tperp = fmin(tperp, stack.z[np] - geometry.zbounds[irz]);
    }
    
    return tperp;
}
/******************************************************************************/

/******************************************************************************/
/* Source definitions */
const int MXEBIN = 200;     // number of energy bins of spectrum
const int INVDIM = 1000;    // number of bins in inverse CDF

struct Source {
    int nmed;                   // number of media in phantom file
    int spectrum;               // 0 : monoenergetic, 1 : spectrum
    int charge;                 // 0 : photons, -1 : electron, +1 : positron
    
    /* For monoenergetic source */
    double energy;
    
    /* For spectrum */
    double deltak;              // number of elements in inverse CDF
    double *cdfinv1;            // energy value of bin
    double *cdfinv2;            // prob. that particle has energy xi
    
    /* Beamlets shape information */
    int nbeamlets;               // number of beamlets per beam
    int *ibeam;                  // index of beam per beamlet
    
    double *xsource;           // coordinates of the source of each beam
    double *ysource;          
    double *zsource;          
        
    double *xcorner;           // coordinates of the bixel corner
    double *ycorner;           
    double *zcorner;  
    
    double *xside1;           // coordinates of the first side of bixel
    double *yside1;           
    double *zside1;
    
    double *xside2;           // coordinates of the second side of bixel
    double *yside2;           
    double *zside2;
        
};
struct Source source;

void initSource() {
    
    /* Get spectrum file path from input data */
    char spectrum_file[128];
    char buffer[BUFFER_SIZE];
    
    source.spectrum = 1;    /* energy spectrum as default case */
    
    /* First check of spectrum file was given as an input */
    if (getInputValue(buffer, "spectrum file") != 1) {
        mexPrintf("Can not find 'spectrum file' key on input file.\n");
        mexPrintf("Switch to monoenergetic case.\n");
        source.spectrum = 0;    /* monoenergetic source */
    }
    
    if (source.spectrum) {
        removeSpaces(spectrum_file, buffer);
        
        /* Open .source file */
        FILE *fp;
        
        if ((fp = fopen(spectrum_file, "r")) == NULL) {
            mexPrintf("Unable to open file: %s\n", spectrum_file);
            exit(EXIT_FAILURE);
        }
        
        mexPrintf("Path to spectrum file : %s\n", spectrum_file);
        
        /* Read spectrum file title */
        fgets(buffer, BUFFER_SIZE, fp);
        mexPrintf("Spectrum file title: %s", buffer);
        
        /* Read number of bins and spectrum type */
        double enmin;   /* lower energy of first bin */
        int nensrc;     /* number of energy bins in spectrum histogram */
        int imode;      /* 0 : histogram counts/bin, 1 : counts/MeV*/
        
        fgets(buffer, BUFFER_SIZE, fp);
        sscanf(buffer, "%d %lf %d", &nensrc, &enmin, &imode);
        
        if (nensrc > MXEBIN) {
            mexPrintf("Number of energy bins = %d is greater than max allowed = "
                   "%d. Increase MXEBIN macro!\n", nensrc, MXEBIN);
            exit(EXIT_FAILURE);
        }
        
        /* upper energy of bin i in MeV */
        double *ensrcd = malloc(nensrc*sizeof(double));
        /* prob. of finding a particle in bin i */
        double *srcpdf = malloc(nensrc*sizeof(double));
        
        /* Read spectrum information */
        for (int i=0; i<nensrc; i++) {
            fgets(buffer, BUFFER_SIZE, fp);
            sscanf(buffer, "%lf %lf", &ensrcd[i], &srcpdf[i]);
        }
        mexPrintf("Have read %d input energy bins from spectrum file.\n", nensrc);
        
        if (imode == 0) {
            mexPrintf("Counts/bin assumed.\n");
        }
        else if (imode == 1) {
            mexPrintf("Counts/MeV assumed.\n");
            srcpdf[0] *= (ensrcd[0] - enmin);
            for(int i=1; i<nensrc; i++) {
                srcpdf[i] *= (ensrcd[i] - ensrcd[i - 1]);
            }
        }
        else {
            mexPrintf("Invalid mode number in spectrum file.");
            exit(EXIT_FAILURE);
        }
        
        double ein = ensrcd[nensrc - 1];
        mexPrintf("Energy ranges from %f to %f MeV\n", enmin, ein);
        
        /* Initialization routine to calculate the inverse of the
         cumulative probability distribution that is used during execution to
         sample the incident particle energy. */
        double *srccdf = malloc(nensrc*sizeof(double));
        
        srccdf[0] = srcpdf[0];
        for (int i=1; i<nensrc; i++) {
            srccdf[i] = srccdf[i-1] + srcpdf[i];
        }
        
        double fnorm = 1.0/srccdf[nensrc - 1];
        double binsok = 0.0;
        source.deltak = INVDIM; /* number of elements in inverse CDF */
        double gridsz = 1.0f/source.deltak;
        
        for (int i=0; i<nensrc; i++) {
            srccdf[i] *= fnorm;
            if (i == 0) {
                if (srccdf[0] <= 3.0*gridsz) {
                    binsok = 1.0;
                }
            }
            else {
                if ((srccdf[i] - srccdf[i - 1]) < 3.0*gridsz) {
                    binsok = 1.0;
                }
            }
        }
        
        if (binsok != 0.0) {
            mexPrintf("Warning!, some of normalized bin probabilities are "
                   "so small that bins may be missed.\n");
        }

        /* Calculate cdfinv. This array allows the rapid sampling for the
         energy by precomputing the results for a fine grid. */
        source.cdfinv1 = malloc(source.deltak*sizeof(double));
        source.cdfinv2 = malloc(source.deltak*sizeof(double));
        double ak;
        
        for (int k=0; k<source.deltak; k++) {
            ak = (double)k*gridsz;
            int i;
            
            for (i=0; i<nensrc; i++) {
                if (ak <= srccdf[i]) {
                    break;
                }
            }
            
            /* We should fall here only through the above break sentence. */
            if (i != 0) {
                source.cdfinv1[k] = ensrcd[i - 1];
            }
            else {
                source.cdfinv1[k] = enmin;
            }
            source.cdfinv2[k] = ensrcd[i] - source.cdfinv1[k];
            
        }
        
        /* Cleaning */
        fclose(fp);
        free(ensrcd);
        free(srcpdf);
        free(srccdf);
    }
    else {  /* monoenergetic source */
        if (getInputValue(buffer, "mono energy") != 1) {
            mexPrintf("Can not find 'mono energy' key on input file.\n");
            exit(EXIT_FAILURE);
        }
        source.energy = atof(buffer);
        mexPrintf("%f monoenergetic source\n", source.energy);
        
    }
    
    /* Parse data of the beamlets */
    unsigned int nfields;
    mxArray *tmp_fieldpointer;

    tmp_fieldpointer = mxGetField(mcSrc,0,"nBixels");
    nfields = mxGetScalar(tmp_fieldpointer);
    source.nbeamlets = nfields;
    
    mexPrintf("%s%d\n", "Total Number of Beamlets:", source.nbeamlets);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"iBeam");
    const double* iBeamPerBeamlet = mxGetPr(tmp_fieldpointer);
    
    source.ibeam = (int*) malloc(source.nbeamlets*sizeof(int));
    for(int i=0; i<source.nbeamlets; i++) {
        source.ibeam[i] = (int) iBeamPerBeamlet[i] - 1; // C indexing style
    }
        
    tmp_fieldpointer = mxGetField(mcSrc,0,"xSource");
    source.xsource = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"ySource");
    source.ysource = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"zSource");
    source.zsource = mxGetPr(tmp_fieldpointer);
            
    tmp_fieldpointer = mxGetField(mcSrc,0,"xCorner");
    source.xcorner = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"yCorner");
    source.ycorner = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"zCorner");
    source.zcorner = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"xSide1");
    source.xside1 = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"ySide1");
    source.yside1 = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"zSide1");
    source.zside1 = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"xSide2");
    source.xside2 = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"ySide2");
    source.yside2 = mxGetPr(tmp_fieldpointer);
    
    tmp_fieldpointer = mxGetField(mcSrc,0,"zSide2");
    source.zside2 = mxGetPr(tmp_fieldpointer);    
    
    return;
}

void cleanSource() {
    
    /* Memory related to the beamlets is freed within Matlab */
    free(source.cdfinv1);
    free(source.cdfinv2);
    
    return;
}

/******************************************************************************/
/* Scoring definitions */
struct Score {
    double ensrc;               // total energy from source
    double *endep;              // 3D dep. energy matrix per batch
    
    /* The following variables are needed for statistical analysis. Their
     values are accumulated across the simulation */
    double *accum_endep;        // 3D deposited energy matrix
    double *accum_endep2;       // 3D square deposited energy
};
struct Score score;

void initScore() {
    
    int gridsize = geometry.isize*geometry.jsize*geometry.ksize;
    
    score.ensrc = 0.0;
    
    /* Region with index 0 corresponds to region outside phantom */
    score.endep = malloc((gridsize + 1)*sizeof(double));
    score.accum_endep = malloc((gridsize + 1)*sizeof(double));
    score.accum_endep2 = malloc((gridsize + 1)*sizeof(double));
    
    /* Initialize all arrays to zero */
    memset(score.endep, 0.0, (gridsize + 1)*sizeof(double));
    memset(score.accum_endep, 0.0, (gridsize + 1)*sizeof(double));
    memset(score.accum_endep2, 0.0, (gridsize + 1)*sizeof(double));
    
    return;
}

void cleanScore() {
    
    free(score.endep);
    free(score.accum_endep);
    free(score.accum_endep2);
    
    return;
}

void ausgab(double edep) {
    
    int np = stack.np;
    int irl = stack.ir[np];
    double endep = stack.wt[np]*edep;
        
    /* Deposit particle energy on spot */
    #pragma omp atomic
    score.endep[irl] += endep;
    
    return;
}

void accumEndep() {
    
    int gridsize = geometry.isize*geometry.jsize*geometry.ksize;
    
    /* Accumulate endep and endep squared for statistical analysis */
    double edep = 0.0;
    
    int irl = 0;
    
    #pragma omp parallel for firstprivate(edep)
    for (irl=0; irl<gridsize + 1; irl++) {
        edep = score.endep[irl];
        
        score.accum_endep[irl] += edep;
        score.accum_endep2[irl] += pow(edep, 2.0);
    }
    
    /* Clean scoring array */
    memset(score.endep, 0.0, (gridsize + 1)*sizeof(double));
    
    return;
}

void accumulateResults(int iout, int nhist, int nbatch)
{
    int irl;
    int imax = geometry.isize;
    int ijmax = geometry.isize*geometry.jsize;
    double endep, endep2, unc_endep;

    /* Calculate incident fluence */
    double inc_fluence = (double)nhist;
    double mass;
    int iz;

    #pragma omp parallel for private(irl,endep,endep2,unc_endep,mass)
    for (iz=0; iz<geometry.ksize; iz++) {
        for (int iy=0; iy<geometry.jsize; iy++) {
            for (int ix=0; ix<geometry.isize; ix++) {
                irl = 1 + ix + iy*imax + iz*ijmax;
                endep = score.accum_endep[irl];
                endep2 = score.accum_endep2[irl];
                
                /* First calculate mean deposited energy across batches and its
                 uncertainty */
                endep /= (double)nbatch;
                endep2 /= (double)nbatch;
                
                /* Batch approach uncertainty calculation */
                if (endep != 0.0) {
                    unc_endep = endep2 - pow(endep, 2.0);
                    unc_endep /= (double)(nbatch - 1);
                    
                    /* Relative uncertainty */
                    unc_endep = sqrt(unc_endep)/endep;
                }
                else {
                    endep = 0.0;
                    unc_endep = 0.9999999;
                }
                
                /* We separate de calculation of dose, to give the user the
                 option to output mean energy (iout=0) or deposited dose
                 (iout=1) per incident fluence */
                
                if (iout) {
                    
                    /* Convert deposited energy to dose */
                    mass = (geometry.xbounds[ix+1] - geometry.xbounds[ix])*
                        (geometry.ybounds[iy+1] - geometry.ybounds[iy])*
                        (geometry.zbounds[iz+1] - geometry.zbounds[iz]);
                    
                    /* Transform deposited energy to Gy */
                    mass *= geometry.med_densities[irl-1];
                    endep *= 1.602E-10/(mass*inc_fluence);
                    
                } else {    /* Output mean deposited energy */
                    endep /= inc_fluence;
                }
                
                /* Store output quantities */
                score.accum_endep[irl] = endep;
                score.accum_endep2[irl] = unc_endep;
            }
        }
    }
    
    /* Zero dose in air */
    #pragma omp parallel for private(irl)
    for (iz=0; iz<geometry.ksize; iz++) {
        for (int iy=0; iy<geometry.jsize; iy++) {
            for (int ix=0; ix<geometry.isize; ix++) {
                irl = 1 + ix + iy*imax + iz*ijmax;
                
                if(geometry.med_densities[irl-1] < 0.044) {
                    score.accum_endep[irl] = 0.0;
                    score.accum_endep2[irl] = 0.9999999;
                }
            }
        }
    }
    
    return;
}

void outputResults(char *output_file, int iout, int nhist, int nbatch) {
    
    /* Accumulate the results */
    accumulateResults(iout, nhist,nbatch);
    
    int irl;
    int imax = geometry.isize;
    int ijmax = geometry.isize*geometry.jsize;
    
    /* Output to file */
    char extension[15];
    if (iout) {
        strcpy(extension, ".3ddose");
    } else {
        strcpy(extension, ".3denergy");
    }
    
    /* Get file path from input data */
    char output_folder[128];
    char buffer[BUFFER_SIZE];
    
    if (getInputValue(buffer, "output folder") != 1) {
        mexPrintf("Can not find 'output folder' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    removeSpaces(output_folder, buffer);
    
    /* Make space for the new string */
    char* file_name = malloc(strlen(output_folder) + strlen(output_file) + 
        strlen(extension) + 1);
    strcpy(file_name, output_folder);
    strcat(file_name, output_file); /* add the file name */
    strcat(file_name, extension); /* add the extension */
    
    FILE *fp;
    if ((fp = fopen(file_name, "w")) == NULL) {
        mexPrintf("Unable to open file: %s\n", file_name);
        exit(EXIT_FAILURE);
    }
    
    /* Grid dimensions */
    fprintf(fp, "%5d%5d%5d\n",
            geometry.isize, geometry.jsize, geometry.ksize);
    
    /* Boundaries in x-, y- and z-directions */
    for (int ix = 0; ix<=geometry.isize; ix++) {
        fprintf(fp, "%f ", geometry.xbounds[ix]);
    }
    fprintf(fp, "\n");
    for (int iy = 0; iy<=geometry.jsize; iy++) {
        fprintf(fp, "%f ", geometry.ybounds[iy]);
    }
    fprintf(fp, "\n");
    for (int iz = 0; iz<=geometry.ksize; iz++) {
        fprintf(fp, "%f ", geometry.zbounds[iz]);
    }
    fprintf(fp, "\n");
    
    /* Dose or energy array */
    for (int iz=0; iz<geometry.ksize; iz++) {
        for (int iy=0; iy<geometry.jsize; iy++) {
            for (int ix=0; ix<geometry.isize; ix++) {
                irl = 1 + ix + iy*imax + iz*ijmax;
                fprintf(fp, "%e ", score.accum_endep[irl]);
            }
        }
    }
    fprintf(fp, "\n");
    
    /* Uncertainty array */
    for (int iz=0; iz<geometry.ksize; iz++) {
        for (int iy=0; iy<geometry.jsize; iy++) {
            for (int ix=0; ix<geometry.isize; ix++) {
                irl = 1 + ix + iy*imax + iz*ijmax;
                fprintf(fp, "%f ", score.accum_endep2[irl]);
            }
        }
    }
    fprintf(fp, "\n");
    
    /* Cleaning */
    fclose(fp);
    free(file_name);

    return;
}

/******************************************************************************/
/* Region-by-region definitions */
void initRegions() {
    
    /* +1 : consider region surrounding phantom */
    int nreg = geometry.isize*geometry.jsize*geometry.ksize + 1;
    
    /* Allocate memory for region data */
    region.med = malloc(nreg*sizeof(int));
    region.rhof = malloc(nreg*sizeof(double));
    region.pcut = malloc(nreg*sizeof(double));
    region.ecut = malloc(nreg*sizeof(double));
    
    /* First get global energy cutoff parameters */
    char buffer[BUFFER_SIZE];
    if (getInputValue(buffer, "global ecut") != 1) {
        mexPrintf("Can not find 'global ecut' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    double ecut = atof(buffer);
    
    if (getInputValue(buffer, "global pcut") != 1) {
        mexPrintf("Can not find 'global pcut' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    double pcut = atof(buffer);
    
    /* Initialize transport parameters on each region. Region 0 is outside the
     geometry */
    region.med[0] = VACUUM;
    region.rhof[0] = 0.0;
    region.pcut[0] = 0.0;
    region.ecut[0] = 0.0;
    
    for (int i=1; i<nreg; i++) {
        
        /* -1 : EGS counts media from 1. Substract 1 to get medium index */
        int imed = geometry.med_indices[i - 1] - 1;
        region.med[i] = imed;
        
        if (imed == VACUUM) {
            region.rhof[0] = 0.0F;
            region.pcut[0] = 0.0F;
            region.ecut[0] = 0.0F;
        }
        else {
            if (geometry.med_densities[i - 1] == 0.0F) {
                region.rhof[i] = 1.0;
            }
            else {
                region.rhof[i] =
                    geometry.med_densities[i - 1]/pegs_data.rho[imed];
            }
            
            /* Check if global cut-off values are within PEGS data */
            if (pegs_data.ap[imed] <= pcut) {
                region.pcut[i] = pcut;
            } else {
                mexPrintf("Warning!, global pcut value is below PEGS's pcut value "
                       "%f for medium %d, using PEGS value.\n",
                       pegs_data.ap[imed], imed);
                region.pcut[i] = pegs_data.ap[imed];
            }
            if (pegs_data.ae[imed] <= ecut) {
                region.ecut[i] = ecut;
            } else {
                mexPrintf("Warning!, global pcut value is below PEGS's ecut value "
                       "%f for medium %d, using PEGS value.\n",
                       pegs_data.ae[imed], imed);
            }
        }
    }

    return;
}

void initHistory(int ibeamlet) {
    
    double rnno1;
    double rnno2;
    
    int ijmax = geometry.isize*geometry.jsize;
    int imax = geometry.isize;
    
    /* Initialize first particle of the stack from source data */
    stack.np = 0;
    stack.iq[stack.np] = source.charge;
    
    /* Get primary particle energy */
    double ein = 0.0;
    if (source.spectrum) {
        /* Sample initial energy from spectrum data */
        rnno1 = setRandom();
        rnno2 = setRandom();
        
        /* Sample bin number in order to select particle energy */
        int k = (int)fmin(source.deltak*rnno1, source.deltak - 1.0);
        ein = source.cdfinv1[k] + rnno2*source.cdfinv2[k];
    }
    else {
        /* Monoenergetic source */
        ein = source.energy;
    }
    
    /* Check if the particle is an electron, in such a case add electron
     rest mass energy */
    if (stack.iq[stack.np] != 0) {
        /* Electron or positron */
        stack.e[stack.np] = ein + RM;
    }
    else {
        /* Photon */
        stack.e[stack.np] = ein;
    }
    
    /* Accumulate sampled kinetic energy for fraction of deposited energy
     calculations */
    score.ensrc += ein;
    
    /* Set particle position. First obtain a random position in the rectangle
     defined by the bixel at isocenter*/    
    double xiso = 0.0; 
    double yiso = 0.0;
    double ziso = 0.0;
    
    rnno1 = setRandom();
    rnno2 = setRandom();
    
    xiso = rnno1*source.xside1[ibeamlet] + rnno2*source.xside2[ibeamlet] + 
            source.xcorner[ibeamlet];
    yiso = rnno1*source.yside1[ibeamlet] + rnno2*source.yside2[ibeamlet] + 
            source.ycorner[ibeamlet];
    ziso = rnno1*source.zside1[ibeamlet] + rnno2*source.zside2[ibeamlet] + 
            source.zcorner[ibeamlet];
    
    /* Norm of the resulting vector from the source of current beam to the 
     position of the particle on bixel */
    int ibeam = source.ibeam[ibeamlet];
    double vnorm = sqrt(pow(xiso - source.xsource[ibeam], 2.0) + 
            pow(yiso - source.ysource[ibeam], 2.0) + 
            pow(ziso - source.zsource[ibeam], 2.0));
        
    /* Direction of the particle from position on bixel to beam source*/
    double u = -(xiso - source.xsource[ibeam])/vnorm;
    double v = -(yiso - source.ysource[ibeam])/vnorm;
    double w = -(ziso - source.zsource[ibeam])/vnorm;
    
    /* Calculate the minimum distance from particle position on bixel to 
     phantom boundaries */
    double ustep = 1.0E5;
    double dist;
    
    if(u > 0.0) {
        dist = (geometry.xbounds[geometry.isize]-xiso)/u;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    if(u < 0.0) {
        dist = -(xiso-geometry.xbounds[0])/u;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    
    if(v > 0.0) {
        dist = (geometry.ybounds[geometry.jsize]-yiso)/v;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    if(v < 0.0) {
        dist = -(yiso-geometry.ybounds[0])/v;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    
    if(w > 0.0) {
        dist = (geometry.zbounds[geometry.ksize]-ziso)/w;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    if(w < 0.0) {
        dist = -(ziso-geometry.zbounds[0])/w;
        if(dist < ustep) {
            ustep = dist;
        }        
    }
    
    /* Transport particle from bixel to surface. Adjust particle direction 
     to be incident to phantom surface */
    stack.x[stack.np] = xiso + ustep*u;
    stack.y[stack.np] = yiso + ustep*v;
    stack.z[stack.np] = ziso + ustep*w;
    
    stack.u[stack.np] = -u;
    stack.v[stack.np] = -v;
    stack.w[stack.np] = -w;

    /* For numerical stability, make sure that points are really inside the phantom */
    if(stack.x[stack.np] < geometry.xbounds[0]) {
        stack.x[stack.np] = geometry.xbounds[0] + 2.0*DBL_MIN;
    }
    if(stack.x[stack.np] > geometry.xbounds[geometry.isize]) {
        stack.x[stack.np] = geometry.xbounds[geometry.isize] - 2.0*DBL_MIN;
    }

    if(stack.y[stack.np] < geometry.ybounds[0]) {
        stack.y[stack.np] = geometry.ybounds[0] + 2.0*DBL_MIN;
    }
    if(stack.y[stack.np] > geometry.ybounds[geometry.jsize]) {
        stack.y[stack.np] = geometry.ybounds[geometry.jsize] - 2.0*DBL_MIN;
    }

    if(stack.z[stack.np] < geometry.zbounds[0]) {
        stack.z[stack.np] = geometry.ybounds[0] + 2.0*DBL_MIN;
    }
    if(stack.z[stack.np] > geometry.zbounds[geometry.ksize]) {
      stack.z[stack.np] = geometry.zbounds[geometry.ksize] - 2.0*DBL_MIN;
    }
    
    /* Determine region index of source particle */
    int ix = 0;
    while (geometry.xbounds[ix+1] < stack.x[stack.np]) {
        ix++;
    }
    
    int iy = 0;
    while (geometry.ybounds[iy+1] < stack.y[stack.np]) {
        iy++;
    }
    
    int iz = 0;
    while (geometry.zbounds[iz+1] < stack.z[stack.np]) {
        iz++;
    }
    
    stack.ir[stack.np] = 1 + ix + iy*imax + iz*ijmax;
          
    /* Set statistical weight and distance to closest boundary*/
    stack.wt[stack.np] = 1.0;
    stack.dnear[stack.np] = 0.0;
    
    return;
}

/******************************************************************************/
/* omc_matrad main function */
void mexFunction (int nlhs, mxArray *plhs[],    // output of the function
    int nrhs, const mxArray *prhs[])            // input of the function
{
    
    /* Execution time measurement */
    double tbegin;
    tbegin = omc_get_time();
    
    /* Parsing program options */

    /* Check for proper number of input and output arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:invalidNumInputs",
                "Two or three input arguments required.");
    }
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "matRad:matRad_ompInterface:invalidNumOutputs",
                "Too many output arguments.");
    }

    parseInput(nrhs, prhs);

    /* Get information of OpenMP environment */
#ifdef _OPENMP
    int omp_size = omp_get_num_procs();
    mexPrintf("Number of OpenMP threads: %d\n", omp_size);
    omp_set_num_threads(omp_size);
#else
    mexPrintf("ompMC compiled without OpenMP support. Serial execution.\n");
#endif
    
    /* Read geometry information from matRad and initialize geometry */
    initPhantom();
    
    /* With number of media and media names initialize the medium data */
    initMediaData();
    
    /* Initialize radiation source */
    initSource();
    
    /* Initialize data on a region-by-region basis */
    initRegions();
    
    /* Initialize VRT data */
    initVrt();
    
    /* Preparation of scoring struct */
    initScore();

    #pragma omp parallel
    {
      /* Initialize random number generator */
      initRandom();

      /* Initialize particle stack */
      initStack();
    }

    /* Shower call */
    
    /* Get number of histories, statistical batches and splitting factor */
    char buffer[BUFFER_SIZE];
    if (getInputValue(buffer, "ncase") != 1) {
        mexPrintf("Can not find 'ncase' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    int nhist = atoi(buffer);
    
    if (getInputValue(buffer, "nbatch") != 1) {
        mexPrintf("Can not find 'nbatch' key on input file.\n");
        exit(EXIT_FAILURE);
    }
    int nbatch = atoi(buffer); 
    
    if (nhist/nbatch == 0) {
        nhist = nbatch;
    }
    
    int nperbatch = nhist/nbatch;
    nhist = nperbatch*nbatch;
    
    int gridsize = geometry.isize*geometry.jsize*geometry.ksize;
    
    mexPrintf("Total number of particle histories: %d\n", nhist);
    mexPrintf("Number of statistical batches: %d\n", nbatch);
    mexPrintf("Histories per batch: %d\n", nperbatch);

    if (getInputValue(buffer, "relative dose threshold") != 1) {
        mexPrintf("Can not find 'relative dose threshold' key on input file.\n");
        exit(EXIT_FAILURE);
    }    
    double relDoseThreshold = atof(buffer);

    mexPrintf("Using a relative dose cut-off of %f\n",relDoseThreshold);
    
    /* Use Matlab waitbar to show execution progress */
    mxArray* waitbarHandle = 0;                             // the waitbar handle does not exist yet
	mxArray* waitbarProgress = mxCreateDoubleScalar(0.0);   // allocate a double scalar for the progress
	mxArray* waitbarMessage = mxCreateString("calculate dose influence matrix for photons (ompMC) ...");    // allocate a string for the message
	
	mxArray* waitbarInputs[3];  // array of waitbar inputs
    mxArray* waitbarOutput[1];  // pointer to waitbar output

	waitbarInputs[0] = waitbarProgress; 
	waitbarInputs[1] = waitbarMessage;
	
	/* Create the waitbar with h = waitbar(progress,message); */
    int status = mexCallMATLAB(1, waitbarOutput, 2, waitbarInputs, "waitbar");

    waitbarHandle = waitbarOutput[0];

    /* Create output matrix */
    mwSize nCubeElements = geometry.isize*geometry.jsize*geometry.ksize;
    double percentage_steps = 0.01;             // steps in which the sparse matrix is allocated
    double percent_sparse = percentage_steps;   // initial percentage to allocate memory for

    mwSize nzmax = (mwSize) ceil((double)nCubeElements*(double)source.nbeamlets*percent_sparse);
    plhs[0] = mxCreateSparse(nCubeElements,source.nbeamlets,nzmax,mxREAL);

    double *sr  = mxGetPr(plhs[0]);
    mwIndex *irs = mxGetIr(plhs[0]);
    mwIndex *jcs = mxGetJc(plhs[0]);
    mwIndex linIx = 0;
    jcs[0] = 0;
    
    double progress = 0.0;

    /* Execution time up to this point */
    mexPrintf("Execution time up to this point : %8.2f seconds\n",
           (omc_get_time() - tbegin));
    
    for(int ibeamlet=0; ibeamlet<source.nbeamlets; ibeamlet++) {
        for (int ibatch=0; ibatch<nbatch; ibatch++) {            
            int ihist;

            #pragma omp parallel for schedule(dynamic)
            for (ihist=0; ihist<nperbatch; ihist++) {
                /* Initialize particle history */
                initHistory(ibeamlet);
                
                /* Start electromagnetic shower simulation */
                shower();
            }
            
            /* Accumulate results of current batch for statistical analysis */
            accumEndep();

            progress = ((double)ibeamlet + (double)(ibatch+1)/nbatch)/source.nbeamlets;
            (*mxGetPr(waitbarProgress)) = progress;

            if (waitbarOutput && waitbarHandle) {              
              waitbarInputs[0] = waitbarProgress;
              waitbarInputs[1] = waitbarHandle;
              waitbarInputs[2] = waitbarMessage;
              status = mexCallMATLAB(0, waitbarOutput, 2, waitbarInputs, "waitbar");
            }
        }

        /* Output of results for current beamlet */
        int iout = 1;   /* i.e. deposit mean dose per particle fluence */
        accumulateResults(iout, nhist, nbatch);

        /* Get maximum value to apply threshold */
        double doseMax = 0.0;
        for (int irl=1; irl < gridsize+1; irl++) {
            if (score.accum_endep[irl] > doseMax) {
                doseMax = score.accum_endep[irl];
            }
        }
        double thresh = doseMax*relDoseThreshold;
        /* Count values above threshold */
        mwSize nnz = 0; //Number of nonzeros in the dose cube

        for (int irl=1; irl < gridsize+1; irl++) {        
            if (score.accum_endep[irl] > thresh) {
                nnz++;
            }                
        }

        /* Check if we need to reallocate for sparse matrix */
        if ((linIx + nnz) > nzmax) {
            int oldnzmax = nzmax;
            percent_sparse += percentage_steps;
            nzmax = (mwSize) ceil((double)nCubeElements*(double)source.nbeamlets*percent_sparse);
            
            /* Make sure nzmax increases at least by 1. */
            if (oldnzmax == nzmax) {
                nzmax++;
            }                

            /* Check that the new nmax is large enough and if not, also adjust 
            the percentage_steps since we seem to have set it too small for this 
            particular use case */
            if (nzmax < (linIx + nnz)) {
                nzmax = linIx + nnz;
                percent_sparse = (double)nzmax/nCubeElements;
                percentage_steps = percent_sparse;
            }

            if (verbose_flag) {
                mexPrintf("Reallocating Sparse Matrix from nzmax=%d to nzmax=%d\n", oldnzmax, nzmax);
            }                
            
            /* Set new nzmax and reallocate more memory */
            mxSetNzmax(plhs[0], nzmax);
            mxSetPr(plhs[0], (double *) mxRealloc(sr, nzmax*sizeof(double)));
            mxSetIr(plhs[0], (mwIndex *)  mxRealloc(irs, nzmax*sizeof(mwIndex)));
            
            /* Use the new pointers */
            sr  = mxGetPr(plhs[0]);
            irs = mxGetIr(plhs[0]);
        }
        
        for (int irl=1; irl < gridsize+1; irl++) {        
            if (score.accum_endep[irl] > thresh) {            
                sr[linIx] = score.accum_endep[irl];
                irs[linIx] = irl-1;
                linIx++;
            }
        }
        
        jcs[ibeamlet+1] = linIx;
        
        /* Reset accum_endep for following beamlet */
        memset(score.accum_endep, 0.0, (gridsize + 1)*sizeof(double));                
		progress = (double) (ibeamlet+1) / (double) source.nbeamlets;		
        (*mxGetPr(waitbarProgress)) = progress;

		/* Update the waitbar with waitbar(hWaitbar,progress); */
        if (waitbarOutput && waitbarHandle) {
            waitbarInputs[0] = waitbarProgress;
            waitbarInputs[1] = waitbarHandle;		
            waitbarInputs[2] = waitbarMessage;
            status = mexCallMATLAB(0, waitbarOutput, 2, waitbarInputs, "waitbar");
        }
    }

    mxDestroyArray (waitbarProgress);
	mxDestroyArray (waitbarMessage);
    if (waitbarOutput && waitbarHandle) {
        waitbarInputs[0] = waitbarHandle;		
        status = mexCallMATLAB(0,waitbarOutput,1, waitbarInputs,"close") ;
        mxDestroyArray(waitbarHandle);
	}
    
    mexPrintf("Sparse MC Dij has %d (%f percent) elements!\n", linIx, 
        (double)linIx/((double)nCubeElements*(double)source.nbeamlets));
    
    /* Truncate the matrix to the exact size by reallocation */
    mxSetNzmax(plhs[0], linIx);
    mxSetPr(plhs[0], mxRealloc(sr, linIx*sizeof(double)));
    mxSetIr(plhs[0], mxRealloc(irs, linIx*sizeof(int)));
    
    sr  = mxGetPr(plhs[0]);
    irs = mxGetIr(plhs[0]);    
    
    /* Print some output and execution time up to this point */
    mexPrintf("Simulation finished\n");
    mexPrintf("Execution time up to this point : %8.2f seconds\n",
           (omc_get_time() - tbegin));    
       
    /* Cleaning */
    cleanPhantom();
    cleanPhoton();
    cleanRayleigh();
    cleanPair();
    cleanElectron();
    cleanMscat();
    cleanSpin();
    cleanRegions();
    cleanScore();
    cleanSource();

    #pragma omp parallel
    {
      cleanRandom();
      cleanStack();
    }

    /* Get total execution time */
    mexPrintf("Total execution time : %8.5f seconds\n",
           (omc_get_time() - tbegin));
    
    return;
    
}
