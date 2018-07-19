function [thirdMom] = matRad_calcThirdRangeMom(numOfSpots,numComp,w,mMeanA,mMeanB,mMeanC,...
                                                                    mWidthA,mWidthB,mWidthC,...
                                                                    mWeightA,mWeightB,mWeightC,...
                                                                    voxelPos,mCovariance) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the third raw moment. 
% 
% call
% [thirdMom] = matRad_calcThirdRangeMom(numOfSpots,numComp,w,mMeanA,mMeanB,mMeanC,...
%                                                                     mWidthA,mWidthB,mWidthC,...
%                                                                     mWeightA,mWeightB,mWeightC,...
%                                                                     voxelPos,mCovariance)
%
% input
%   numOfSpots           total number of spots
%   numComp              total number of Gaussian components e.g. 10 or 13
%   mMeanA-mMeanC:      four mean vectors of spot components j o q n
%   mWidthA-mWidthC:     four width vectors of spot components j o q n
%   mWeightA-mWeightC:   four weight vectors of spot components j o q n
%   voxelPos:            voxel position 
%   mSysCovRadDept:      full covariance matrix of range error
%
% output
%   third:               third central moment
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
assert(isa(mWidthA,'double'));
assert(isa(mWidthB,'double'));
assert(isa(mWidthC,'double'));
assert(isa(mWeightA,'double'));
assert(isa(mWeightB,'double'));
assert(isa(mWeightC,'double'));
assert(isa(voxelPos,'double'));
assert(isa(mCovariance,'double'));

coder.varsize('numOfSpots',      [1 1],[0 0]);
coder.varsize('numComp',         [1 1],[0 0]);
coder.varsize('w',               [60 1],[1 0]);
coder.varsize('mMeanA',          [15 60],[1 1]);
coder.varsize('mMeanB',          [15 60],[1 1]);
coder.varsize('mMeanC',          [15 60],[1 1]);
coder.varsize('mWidthA',         [15 60],[1 1]);
coder.varsize('mWidthB',         [15 60],[1 1]);
coder.varsize('mWidthC',         [15 60],[1 1]);
coder.varsize('mWeightA',        [15 60],[1 1]);
coder.varsize('mWeightB',        [15 60],[1 1]);
coder.varsize('mWeightC',        [15 60],[1 1]);
coder.varsize('voxelPos',        [1 1],[0 0]);
coder.varsize('mCovariance',  [60 60],[1 1]);


mPSI_jop  = zeros(numOfSpots,numOfSpots,numOfSpots);

for j = 1:numOfSpots
    for o = 1:numOfSpots
        for q = 1:numOfSpots


            A11_vec  = mWidthA(:,j) + mCovariance(j,j);    %       | vA11 sA12 sA13 | vA11 sA12 sA13
            A22_vec  = mWidthB(:,o) + mCovariance(o,o);    %       | sA21 vA22 sA23 | sA21 vA22 sA23
            A33_vec  = mWidthC(:,q) + mCovariance(q,q);    %       | sA31 sA32 vA33 | sA31 sA32 vA33
            A12      = mCovariance(j,o);
            A13      = mCovariance(j,q);
            A23      = mCovariance(o,q);
                                       
            mDet     = reshape(reshape(A11_vec*A22_vec',[],1) * A33_vec',[numComp numComp numComp]) + A12*A23*A13 + A12*A23*A13;
            a        = A13*(ones(numComp,1)); aa = a * A22_vec'; cube1 = reshape(reshape(aa,[],1) * a',[numComp numComp numComp]);
            b        = A23*(ones(numComp,1)); bb = A11_vec * b'; cube2 = reshape(reshape(bb,[],1) * b',[numComp numComp numComp]);
            c        = A12*(ones(numComp,1)); cc = c * c'; cube3 = reshape(reshape(cc,[],1) * A33_vec',[numComp numComp numComp]);
            mDet     = mDet -cube1 - cube2 - cube3;


            sigma   = [mCovariance(j,j)   mCovariance(j,o)  mCovariance(j,q);...
                       mCovariance(o,j)   mCovariance(o,o)  mCovariance(o,q);...
                       mCovariance(q,j)   mCovariance(q,o)  mCovariance(q,q) ];

            dev_j   = voxelPos - mMeanA(:,j);
            dev_o   = voxelPos - mMeanB(:,o);
            dev_q   = voxelPos - mMeanC(:,q);


            tmp = reshape(mWeightA(:,j) * mWeightB(:,o)',[numComp^2 1]);
            vW  = reshape(tmp * mWeightC(:,q)',[numComp numComp numComp]);
            
            %vWref = reshape(reshape(mWeightA(:,j) * mWeightB(:,o)',[numComp^2 1]) * mWeightB(:,q)',[numComp numComp numComp]);
        
            Y = 0; %Yref = 0;
            for J = 1:numComp
                 for O = 1:numComp
                      for Q = 1:numComp  

                        lambda = diag([mWidthA(J,j) mWidthB(O,o) mWidthC(Q,q)]);  

                        LaSi    = (lambda + sigma);
                        InvLaSi = inv(LaSi);

                         Y = Y + (vW(J,O,Q)/((2*pi)^(3/2)*sqrt(mDet(J,O,Q)))) * ...
                             exp(-0.5*([dev_j(J) dev_o(O) dev_q(Q)]) * (InvLaSi) * ([dev_j(J) dev_o(O) dev_q(Q)])');
                          
                          
                         %Yref = Yref + (vWref(J,O,Q) * mvnpdf([dev_j(J) dev_o(O) dev_q(Q)],0,LaSi)); 

%                          mDetRef = [mWidthA(J,j) + mCovariance(j,j)        A12                                     A13;
%                                     A12                                       mWidthB(O,o) + mCovariance(o,o)      A23;
%                                     A13                                       A23                                     mWidthB(Q,q) + mCovariance(q,q)];
%                                 
%                          if ((det(mDetRef)/mDet(J,O,Q)) -1)  >0.01
%                              warning('det are not equal');
%                          end
%                          
%                          if ~isequal(vW(J,O,Q),mWeightA(J,j)  * mWeightB(O,o) * mWeightB(Q,q))
%                              warning('weights are not equal');
%                          end
                      end
                 end
            end

            mPSI_jop(j,o,q) = Y;
            
        end
    end
end

thirdMom = w' *  reshape(reshape(mPSI_jop,[numOfSpots^2 numOfSpots]) * w, [numOfSpots numOfSpots]) * w;    

end
