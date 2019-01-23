function mx = matRad_sparseBeamletsReaderMSsquareXXX(filename,numOfBixels,cubeDim,mask)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsqaure monte carlo photon dose calculation wrapper
%
% call
%   mx = matRad_sparseBeamletsReaderMSsquare(filename,numOfBixels, cubeDim)
%
% input
%   filename:    filename of binary MCsquare dose influence matrix
%   numOfBixels: number of bixels in sparse matrix
%   cubeDim:     dimensions of the dose cube
% output
%   mx:          sparse dose influence matrix in matlab
%
% References
%
%   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX there needs to go a lot of stuff!
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
numOfVoxels = prod(cubeDim);
numOfVoxelsSlice = prod(cubeDim([1 2]));

mx = spalloc(numOfVoxels,numOfBixels,1);

fileHandle = fopen(filename,'r');

for i = 1:numOfBixels
   
    tmpVec = zeros(numOfVoxels,1);
    counter = 0;
    
    nonZeroVoxels = fread(fileHandle,1,'uint32');
    beamID        = fread(fileHandle,1,'uint32');
    layerID       = fread(fileHandle,1,'uint32');
    xCoord        = fread(fileHandle,1,'float');
    yCoord        = fread(fileHandle,1,'float');
    
    while counter < nonZeroVoxels
       
        numOfContValues = fread(fileHandle,1,'uint32');
        
        counter = counter + numOfContValues;
        
        firstIndex = fread(fileHandle,1,'uint32');
        
        values  = fread(fileHandle,numOfContValues,'float');
        ixMCsquare = firstIndex + [1:numOfContValues] + 1;
        
        kSub = ceil(ixMCsquare / numOfVoxelsSlice);
        jSub = ceil( (ixMCsquare - (kSub-1)*numOfVoxelsSlice)/cubeDim(2) );
        iSub = ixMCsquare - (jSub-1)*cubeDim(1) - (kSub-1)*numOfVoxelsSlice;
                    
        % flip image
        iSub = cubeDim(2) - iSub + 1;

        % compute new index
        ixMatRad = jSub + (iSub-1)*cubeDim(2) + (kSub-1)*numOfVoxelsSlice;
        
        % mask
        currMask = mask(ixMatRad);
        
        % fill
        tmpVec(ixMatRad(currMask)) = values(currMask);
        
    end
    
    mx(:,i) = sparse(tmpVec);
    
end

