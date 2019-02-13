/*
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% The original version of this file is part of openREGGUI and can be found 
% at https://openreggui.org/git/open/REGGUI/blob/master/plugins/openMIROpt/
% functions/io/mexSparseBeamletsReader_double_OptROIvoxels.c
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
compile with matlab: mex -largeArrayDims mexSparseBeamletsReader_double_OptROIvoxels.c
run with matlab: Beamlets = mexSparseBeamletsReader('Sparse_Dose.bin', [256 256 71], 8372);
*/

#include "mex.h"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <array>
#include <vector>
#include <cmath>

//typedef std::vector<size_t> ixVec_t;

//Function to reorder the vectors
void reorder_entries(mwIndex* nonZeroIx, double* nonZeroVals, size_t numValues)
{
  //Create a vector for sorting
  std::vector<size_t> sortIx(numValues);
  std::iota(sortIx.begin(), sortIx.end(), 0); //Fills with increasing values
  
  //determine the order by sorting with respect to the voxel indices
  std::sort(sortIx.begin(), sortIx.end(), [nonZeroIx](mwSize i, mwSize j) -> bool { return nonZeroIx[i] < nonZeroIx[j]; });
  
  //Reorder both vectors according to sortIx
  for (int i = 0; i < numValues - 1; ++i)
  {
    //we need to make sure that we put the element in place without overwriting 
    while (i != sortIx[i])
    {
      int alt = sortIx[i];
      std::swap(nonZeroIx[i], nonZeroIx[alt]);
      std::swap(nonZeroVals[i], nonZeroVals[alt]);
      std::swap(sortIx[i], sortIx[alt]);
    }
  }
}


//#include <stdlib.h>
//#include <stdint.h>
//#include <math.h>

/*                  NbrOutputs  Outputs             NbrInputs   Inputs  */
void mexFunction(   int nlhs,   mxArray *plhs[],    int nrhs,   const mxArray *prhs[]){
    
    if(nrhs != 4) 
		mexErrMsgIdAndTxt("matRad_sparseBeamletsReaderMCsquare:nrhs", "4 arguments required, Call matRad_sparseBeamletsReaderMCsquare(filename, doseGridSize, nSpots, voxelIx)");
	
    if(nlhs > 1) 
		mexErrMsgIdAndTxt("matRad_sparseBeamletsReaderMCsquare:nlhs", "Too many output arguments.");
    
    //Check Filename
    if(!mxIsChar(prhs[0])) 
		mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:inputNotString", "filename must be a string/char array");
	
    char *filename = mxArrayToString(prhs[0]);
		
  
    //Check dose grid
    if(mxGetNumberOfElements(prhs[1]) != 3) 
		mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:invalidInput", "Dose grid must be 3D");
    
	  std::array<uint32_t,3> sizeDoseGrid;    
    sizeDoseGrid[0] = (uint32_t) *(mxGetPr(prhs[1])+0);
    sizeDoseGrid[1] = (uint32_t) *(mxGetPr(prhs[1])+1);
    sizeDoseGrid[2] = (uint32_t) *(mxGetPr(prhs[1])+2);
	
	  uint32_t numVoxels = std::accumulate(sizeDoseGrid.begin(), sizeDoseGrid.end(), 1, std::multiplies<uint32_t>());
    //mexPrintf("Number of voxels: %d", numVoxels);
	  uint32_t numVoxelsSlice = sizeDoseGrid[0]*sizeDoseGrid[1];


    //Check numSpots
    if(!mxIsScalar(prhs[2])) 
		mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:invalidInput", "number of spots must be a scalar");
    uint32_t numSpots = (uint32_t) mxGetScalar(prhs[2]);
    
    //Check (logical) voxel indices
    if (!mxIsLogical(prhs[3]) || mxGetN(prhs[3]) != 1) {
        mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:invalidInput", "voxel indicies must be a logical index vector");
    } 
    mxLogical *x = mxGetLogicals(prhs[3]); 
    
    if(mxGetM(prhs[3]) != numVoxels) {
        mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:invalidInput", "size of dose cube and voxel index vector are not the same!");
    }
    
    
    //Now that everything is okay, try to open the file
	  std::ifstream file(filename,std::ios::binary);
    if(!file.is_open()) 
		mexErrMsgIdAndTxt( "matRad_sparseBeamletsReaderMCsquare:fileOpenFailed", "File could not be opened");
	
	
    double percent_sparse = 0.005;
    mwSize nnzMaxEstimate = (mwSize) ((double) numVoxels * (double) numSpots * percent_sparse);
    mwSize oldNnzMaxEstimate;
	
    plhs[0] = mxCreateSparse(numVoxels,numSpots,nnzMaxEstimate,mxREAL);
    double *vals = mxGetPr(plhs[0]); //Value Pointer
    mwIndex *ix = mxGetIr(plhs[0]); //Index Pointer
    mwIndex *jc = mxGetJc(plhs[0]); //Array to keep track of columns
    mwIndex currentNnz = 0;
    
    //For Beamlet information
    uint32_t nBeamletVoxels, iBeam, iLayer; 
    float spotX, spotY;

    //For reading
    uint32_t numCurrentReadVals, numReadValsTot, ixFirst;

    //Data array
    float *data = (float*) mxCalloc(numVoxels, sizeof(float));
    
    //Matrad indexing
    uint32_t ixMatRad, ixMC2, iSub, jSub, kSub;
    
    for(uint32_t runSpot=0; runSpot < numSpots; ++runSpot) 
	  {   
        //Read Beamlet Information
        file.read((char *) &nBeamletVoxels, sizeof(uint32_t));
        file.read((char *) &iBeam, sizeof(uint32_t));
        file.read((char *) &iLayer, sizeof(uint32_t));
        file.read((char *) &spotX, sizeof(float));
        file.read((char *) &spotY, sizeof(float));
    
        numReadValsTot = 0;
        jc[runSpot] = currentNnz;

        while(true) {
            file.read((char *) &numCurrentReadVals, sizeof(uint32_t));
            
            file.read((char *) &ixFirst, sizeof(uint32_t));
            
            file.read((char *) data, numCurrentReadVals * sizeof(float));
            
            for(uint32_t readIx = 0; readIx < numCurrentReadVals; ++readIx){
                
                //If necessary reallocate space for non-zeros
                if(currentNnz >= nnzMaxEstimate){
                    oldNnzMaxEstimate = nnzMaxEstimate;
                    percent_sparse += 0.001;
                    nnzMaxEstimate = (mwSize) std::ceil((double)numVoxels * (double)numSpots * percent_sparse);
                    
                    if (nnzMaxEstimate <= oldNnzMaxEstimate) nnzMaxEstimate = oldNnzMaxEstimate + 1;
                    
                    mxSetNzmax(plhs[0], nnzMaxEstimate);
                    mxSetPr(plhs[0], (double *) mxRealloc(vals, nnzMaxEstimate*sizeof(double)));
                    mxSetIr(plhs[0], (mwIndex *) mxRealloc(ix, nnzMaxEstimate*sizeof(mwIndex)));
                    
                    vals  = mxGetPr(plhs[0]);
                    ix = mxGetIr(plhs[0]);
                }
                
                // get index
                ixMC2 = ixFirst + readIx;
                    
                // get subscripts
                kSub = ixMC2/numVoxelsSlice;
                jSub = (ixMC2 - kSub*numVoxelsSlice)/sizeDoseGrid[1];
                iSub = ixMC2 - jSub*sizeDoseGrid[1] - kSub*numVoxelsSlice;
                    
                // flip image
                iSub = sizeDoseGrid[0] - iSub - 1;

                // compute new index
                ixMatRad = jSub + iSub*sizeDoseGrid[0] + kSub*numVoxelsSlice;
                
                if (x[ixMatRad]) {
                    vals[currentNnz] = (double) data[readIx];
                    ix[currentNnz]   = ixMatRad;
                    currentNnz ++;
                }
            }
            

            numReadValsTot += numCurrentReadVals;
            if(numReadValsTot >= nBeamletVoxels) 
                break;
        }
        uint32_t addedVals = currentNnz - jc[runSpot]; // Number of values actually added
        //This is were the index and value lists for the spot/column starts
        mwIndex* ixBegin = &ix[jc[runSpot]];
        double* valsBegin = &vals[jc[runSpot]];
        
        //Reorders the entries such that they are increasing in the matRAd voxel index
        reorder_entries(ixBegin, valsBegin, addedVals);
    }
    
    jc[numSpots] = currentNnz;
    
    mxFree(data);
    file.close();

    //shrink memory to fit exactly the number of elements
    mxSetNzmax(plhs[0], currentNnz);
    mxSetPr(plhs[0], (double *)  mxRealloc(vals, currentNnz * sizeof(double)));
    mxSetIr(plhs[0], (mwIndex *) mxRealloc(ix  , currentNnz * sizeof(mwIndex)));
}



