function [fourthMom] = matRad_calcFourthRangeMom(numOfSpots,numComp,w,mMeanA,mMeanB,mMeanC,mMeanD,...
                                                                      mWidthA,mWidthB,mWidthC,mWidthD,...
                                                                      mWeightA,mWeightB,mWeightC,mWeightD,...
                                                                      voxelPos,mSysCovRadDept) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the fourth raw moment. 
% 
% call
%   [fourthMom] = matRad_calcFourthRangeMom(numOfSpots,numComp,w,mMeanA,mMeanB,mMeanC,mMeanD, ...
%                                                                mWidthA,mWidthB,mWidthC,mWidthD,...
%                                                                mWeightA,mWeightB,mWeightC,mWeightD,...
%                                                                voxelPos,mSysCovRadDept)
%
% input
%   numOfSpots           total number of spots
%   numComp              total number of Gaussian components e.g. 10 or 13
%   mMeanA-mMeanD:      four mean vectors of spot components j o q n
%   mWidthA-mWidthD:     four width vectors of spot components j o q n
%   mWeightA-mWeightD:   four weight vectors of spot components j o q n
%   voxelPos:            voxel position 
%   mSysCovRadDept:      full covariance matrix of range error
%
% output
%   fourthMom:           fourth central moment
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Hans-Peter Wieser
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(isa(numOfSpots,'double'));
assert(isa(numComp,'double'));
assert(isa(w,'double'));
assert(isa(mMeanA,'double'));
assert(isa(mMeanB,'double'));
assert(isa(mMeanC,'double'));
assert(isa(mMeanD,'double'));
assert(isa(mWidthA,'double'));
assert(isa(mWidthB,'double'));
assert(isa(mWidthC,'double'));
assert(isa(mWidthD,'double'));
assert(isa(mWeightA,'double'));
assert(isa(mWeightB,'double'));
assert(isa(mWeightC,'double'));
assert(isa(mWeightD,'double'));
assert(isa(voxelPos,'double'));
assert(isa(mSysCovRadDept,'double'));

coder.varsize('numOfSpots',      [1 1], [0 0]);
coder.varsize('numComp',         [1 1], [0 0]);
coder.varsize('w',               [60 1],[1 0]);
coder.varsize('mMeanA',          [15 60],[1 1]);
coder.varsize('mMeanB',          [15 60],[1 1]);
coder.varsize('mMeanC',          [15 60],[1 1]);
coder.varsize('mMeanD',          [15 60],[1 1]);
coder.varsize('mWidthA',         [15 60],[1 1]);
coder.varsize('mWidthB',         [15 60],[1 1]);
coder.varsize('mWidthD',         [15 60],[1 1]);
coder.varsize('mWidthD',         [15 60],[1 1]);
coder.varsize('mWeightA',        [15 60],[1 1]);
coder.varsize('mWeightB',        [15 60],[1 1]);
coder.varsize('mWeightC',        [15 60],[1 1]);
coder.varsize('mWeightD',        [15 60],[1 1]);
coder.varsize('voxelPos',        [1 1],[0 0]);
coder.varsize('mSysCovRadDept',  [60 60],[1 1]);


mPSI_joqn = zeros(numOfSpots,numOfSpots,numOfSpots,numOfSpots);

for j = 1:numOfSpots
    for o = 1:numOfSpots
        for q = 1:numOfSpots
           for n = 1:numOfSpots

            sigma = [mSysCovRadDept(j,j)  mSysCovRadDept(j,o)    mSysCovRadDept(j,q)   mSysCovRadDept(j,n);...
                     mSysCovRadDept(o,j)  mSysCovRadDept(o,o)    mSysCovRadDept(o,q)   mSysCovRadDept(o,n);...
                     mSysCovRadDept(q,j)  mSysCovRadDept(q,o)    mSysCovRadDept(q,q)   mSysCovRadDept(q,n);...
                     mSysCovRadDept(n,j)  mSysCovRadDept(n,o)    mSysCovRadDept(n,q)   mSysCovRadDept(n,n)];

            dev_j   = voxelPos - mMeanA(:,j);
            dev_o   = voxelPos - mMeanB(:,o);
            dev_q   = voxelPos - mMeanC(:,q);
            dev_n   = voxelPos - mMeanD(:,n);
            
            vW = reshape(reshape(reshape(mWeightA(:,j) * mWeightB(:,o)',[],1) * mWeightC(:,q)',[],1) * mWeightD(:,n)',[numComp,numComp,numComp,numComp]);
            
            Y = 0;
            for J = 1:numComp
                 for O = 1:numComp
                      for Q = 1:numComp  
                          for N = 1:numComp

                            lambda = diag([mWidthA(J,j) mWidthB(O,o) mWidthC(Q,q) mWidthD(N,n)]);  

                            LaSi    = (lambda + sigma);
                            InvLaSi = inv(LaSi);
                        
                             Y = Y + (vW(J,O,Q,N)/((2*pi)^(4/2)*sqrt(det(LaSi)))) * ...
                                 exp(-0.5*([dev_j(J) dev_o(O) dev_q(Q) dev_n(N)]) * (InvLaSi) * ([dev_j(J) dev_o(O) dev_q(Q) dev_n(N)])');

                          end
                      end
                 end
            end
                        
            mPSI_joqn(j,o,q,n) = Y;
            
           end
        end
    end
end


fourthMom =  w'* reshape(reshape(reshape(mPSI_joqn,[numOfSpots^3 numOfSpots]) * w, [numOfSpots^2 numOfSpots]) * w , [numOfSpots numOfSpots]) * w;    


end
