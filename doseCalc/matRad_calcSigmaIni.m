function [sigmaIni] = matRad_calcSigmaIni(baseData,rays,SSD)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates simultaneously the initial sigma of the beam for
% one or more energies
%
% call
%   sigmaIni = matRad_calcSigmaIni(machine.data,stf(i).ray,stf(i).ray(j).SSD);
%
% input
%   baseData:           'machine.data' file
%   rays:            	'stf.ray' file
%   SSD:                source-surface difference
%
% output
%   sigmaIni:        	initial sigma of the ray at certain energy (or
%                       energies). The data is given in 1xP dimensions,
%                       where 'P' represents the number of different
%                       energies
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is not part of the offical matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds the energies involved in the process and their focus index
[energyVec,uniqueEnergyIdx] = unique([rays.energy]);
tempFocus = [rays.focusIx];
focusIxVec = tempFocus(uniqueEnergyIdx);

% finds energy index for all the energies
[~,energyIxVec] = intersect([baseData.energy],energyVec);

if length(energyVec)==1
    sigmaIni = matRad_interp1(baseData(energyIxVec).initFocus.dist(rays.focusIx,:)',baseData(energyIxVec).initFocus.sigma(rays.focusIx,:)',SSD);
    
else
    % finds standard deviation and distance for all the energies
    focusStruct = [baseData.initFocus];
    sigmaVec = [focusStruct(energyIxVec).sigma];
    distVec = [focusStruct(energyIxVec).dist];
    
    dim = size(focusStruct(1).sigma,2);
    
    % repeats every index for a number of times equal to size of sigma
%     repFocusIxVec = kron(focusIxVec, ones(1,dim));
    repFocusIxVec = repelem(focusIxVec,dim);
    
    % Take the values of sigma and distance that we need for calculation]
    % The elements I want are the ones which corrispond to index in the same
    % position in both index vectors, i.e. the diagonal of the resulting matrix
    sigmaVec = sigmaVec(repFocusIxVec,[1:size(sigmaVec,2)]);
    sigmaVec = diag(sigmaVec);
    distVec = distVec(repFocusIxVec,[1:size(distVec,2)]);
    distVec = diag(distVec);
    
    % reshape because I want as output a vector of the correct sigmas and dist,
    % so I interpolate 'dim' elements at the time
    distMat = reshape(distVec,dim,[]);
    sigMat = reshape(sigmaVec,dim,[]);
    
    [x,~]=meshgrid(1:size(distMat,2),1:1);
    [x1,~]=meshgrid(1:size(distMat,2),1:size(distMat,1));
    
    % griddata produces interpolation error in octave
    % sigmaIni = griddata(distMat,x1, sigMat, repmat(SSD,1,size(distMat,2)),x);
    sigmaIni = interp2(x1,distMat, sigMat,x, repmat(SSD,1,size(distMat,2)));
end
