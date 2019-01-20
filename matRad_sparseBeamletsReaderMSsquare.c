/*
compile with matlab: mex -largeArrayDims mexSparseBeamletsReader_double_OptROIvoxels.c
run with matlab: Beamlets = mexSparseBeamletsReader('Sparse_Dose.bin', [256 256 71], 8372);
*/

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/*                  NbrOutputs  Outputs             NbrInputs   Inputs  */
void mexFunction(   int nlhs,   mxArray *plhs[],    int nrhs,   const mxArray *prhs[]){
    
    if(nrhs < 4) mexErrMsgIdAndTxt("mexSparseReader:nrhs", "This function takes 3 input arguments: file name, Dose grid size, Number of spots and voxel indices");
    else if(nlhs > 1) mexErrMsgIdAndTxt("mexSparseReader:nlhs", "Too many output arguments.");
    
    /* first input must be a string */
    if(mxIsChar(prhs[0]) != 1) mexErrMsgIdAndTxt( "mexSparseReader:inputNotString", "First input must be a string.");
    if(mxGetM(prhs[0])!=1) mexErrMsgIdAndTxt( "mexSparseReader:inputNotVector", "First input must be a row vector.");
    char *FileName;
    FileName = mxArrayToString(prhs[0]);
    if(FileName == NULL) mexErrMsgIdAndTxt( "mexSparseReader:InputConversionFailed", "Could not convert input to string.");
  
    /* second input must be a numeric array */
    if(mxGetNumberOfElements(prhs[1]) != 3) mexErrMsgIdAndTxt( "mexSparseReader:InvalidInputArray", "Dose grid sizes must contain 3 elements.");
    double DoseGridSize[3];
    DoseGridSize[0] = *(mxGetPr(prhs[1])+0);
    DoseGridSize[1] = *(mxGetPr(prhs[1])+1);
    DoseGridSize[2] = *(mxGetPr(prhs[1])+2);
    int NbrVoxels = (int)DoseGridSize[0]*DoseGridSize[1]*DoseGridSize[2];

    /* third input must be numeric */
    if(mxIsNumeric(prhs[2]) != 1) mexErrMsgIdAndTxt( "mexSparseReader:InvalidNbrSpots", "Third argument (number of spots) must be numeric.");
    int NbrSpots = (int)mxGetScalar(prhs[2]);
    
        /* fourth input must be a logical vector*/
    if (! mxIsLogical(prhs[3]) || mxGetN(prhs[3]) != 1 ) {
        mexErrMsgTxt("Fourth argument should be an logical vector (double precision).");
    } 
    mxLogical *x;
    x = mxGetLogicals(prhs[3]); /*input vector*/
    
    if(mxGetM(prhs[3]) != NbrVoxels) {
        mexErrMsgTxt("Dimensions of matrix and vector do not match."); /* check that the vector has the size equal to the number of voxels*/
    }
    
    
    /* Create the sparse matrix */
    double percent_sparse = 0.009;
    mwSize NonZeroMax = (mwSize)ceil((double)NbrVoxels * (double)NbrSpots * percent_sparse);
    mwSize oldNonZeroMax;
    plhs[0] = mxCreateSparse(NbrVoxels,NbrSpots,NonZeroMax,mxREAL);
    double *Sparse_data = mxGetPr(plhs[0]);
    mwIndex *Sparse_index = mxGetIr(plhs[0]);
    mwIndex *Sparse_JC = mxGetJc(plhs[0]); /* Sparse_index of the first Non Zero value in each column */
    mwIndex NbrElements = 0;
    
    FILE *fid = NULL;
    fid = fopen(FileName, "rb");
    
    if(fid == NULL) mexErrMsgIdAndTxt( "mexSparseReader:UnableToOpenFile", "Unable to open the binary file.");
    
    uint32_t NonZeroVoxels, BeamID, LayerID, NbrContinuousValues, ReadVoxels, FirstIndex;
    float xcoord, ycoord;
    size_t NbrBytesRead;
    int i, j;
    float *data = mxCalloc(NbrVoxels, sizeof(float));
    
    for(i=0; i<NbrSpots; i++){
        
        /* mexPrintf("Spot %d / %d \n", i+1, NbrSpots); */
        
        NbrBytesRead = fread(&NonZeroVoxels, sizeof(uint32_t), 1, fid);
        NbrBytesRead = fread(&BeamID, sizeof(uint32_t), 1, fid);
        NbrBytesRead = fread(&LayerID, sizeof(uint32_t), 1, fid);
        NbrBytesRead = fread(&xcoord, sizeof(float), 1, fid);
        NbrBytesRead = fread(&ycoord, sizeof(float), 1, fid);
    
        ReadVoxels = 0;
        Sparse_JC[i] = NbrElements;
  
        while(1){
            NbrBytesRead = fread(&NbrContinuousValues, sizeof(uint32_t), 1, fid);
            ReadVoxels += NbrContinuousValues;
            
            NbrBytesRead = fread(&FirstIndex, sizeof(uint32_t), 1, fid);
            FirstIndex = NbrVoxels - FirstIndex;
            
            NbrBytesRead = fread(data, sizeof(float), NbrContinuousValues, fid);
            
            for(j=0; j < NbrContinuousValues; j++){
                
                if(NbrElements >= NonZeroMax){
                    oldNonZeroMax = NonZeroMax;
                    percent_sparse += 0.001;
                    NonZeroMax = (mwSize)ceil((double)NbrVoxels * (double)NbrSpots * percent_sparse);
                    
                    mexPrintf("Realloc %d -> %d \n", oldNonZeroMax, NonZeroMax); 
                    
                    /* make sure nzmax increases atleast by 1 */
                    if (NonZeroMax <= oldNonZeroMax) NonZeroMax = oldNonZeroMax + 1;
                    
                    mxSetNzmax(plhs[0], NonZeroMax);
                    mxSetPr(plhs[0], mxRealloc(Sparse_data, NonZeroMax*sizeof(double)));
                    mxSetIr(plhs[0], mxRealloc(Sparse_index, NonZeroMax*sizeof(mwIndex)));
                    
                    Sparse_data  = mxGetPr(plhs[0]);
                    Sparse_index = mxGetIr(plhs[0]);
                }
                
                if (x[FirstIndex-j-1]){
                    Sparse_data[NbrElements] = (double)data[j];
                    Sparse_index[NbrElements] = FirstIndex-j-1;
                    NbrElements ++;
                }
            }
            
            if(ReadVoxels >= NonZeroVoxels) break;
        }
  
    }
    
    Sparse_JC[NbrSpots] = NbrElements;
    
    mxFree(data);
    fclose(fid);
}
