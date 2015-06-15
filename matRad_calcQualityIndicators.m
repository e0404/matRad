function QI = matRad_calcQualityIndicators(d,cst,refGy,refVol)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dvh calculation
% 
% call
%   matRad_calcDVH(d,cst,lineStyleIndicator)
%
% input
%   d:                  dose cube
%   cst:                matRad cst struct
%   refGy:              array of dose values used for V_XGy calculation
%   refVol:             array of volumes (0-100) used for D_X calculation
%
% output
%   various quality indicators like CI, HI (for targets) and DX, VX within 
%   a structure set   
%
% References
%   van't Riet et. al., IJROBP, 1997 Feb 1;37(3):731-6.
%   Kataria et. al., J Med Phys. 2012 Oct-Dec; 37(4): 207–213.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin < 3)
    refGy = [40 50 60];
    refVol = [2 5 95 98];
end

numOfVois = size(cst,1);

%% calculate QIs per VOI
for i = 1:numOfVois  
    indices     = cst{i,4};
    numOfVoxels = numel(indices);
    
    %% get Dose, dose is sorted to simplify calculations
    if sum(strcmp(fieldnames(d),'RBEWeightedDose')) > 0
        relevantDose = d.RBEWeightedDose;
        doseInVoi   = sort(d.RBEWeightedDose(indices));
    else
        relevantDose = d.physicalDose;
        doseInVoi   = sort(d.physicalDose(indices));
    end
        
    %"easy ones"
    QI(i).mean = mean(doseInVoi);
    QI(i).std = std(doseInVoi);
    QI(i).max = doseInVoi(end);
    QI(i).min = doseInVoi(1);
    
    DX = @(x) doseInVoi(round((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
    
    %create runtime structure fields.. can this be written more elegant?
    for(runDX = 1:numel(refVol))
        QI(i).(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));
    end
    for(runVX = 1:numel(refGy))
        QI(i).(strcat('V',num2str(refGy(runVX)))) = VX(refGy(runVX));
    end
    
    
    if strcmp(cst{i,3},'TARGET') > 0       
        %% Conformity Index
        VTarget95 = numel(doseInVoi(doseInVoi >= 0.95*cst{i,6}.parameter(2)));   %number of target voxels recieving dose >= 0.95 dPres
        VTreated95 = numel(relevantDose(relevantDose >= 0.95*cst{i,6}.parameter(2)));  %number of all voxels recieving dose >= 0.95 dPres ("treated volume")
        QI(i).CI = VTarget95^2/(numOfVoxels * VTreated95);
        
        %% Homogeneity Index (one out of many)        
        QI(i).HI = (DX(5) - DX(95))/cst{i,6}.parameter(2) * 100;
    end
end
