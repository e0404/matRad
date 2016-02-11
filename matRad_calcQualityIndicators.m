function result = matRad_calcQualityIndicators(result,cst,refGy,refVol)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad QI calculation
% 
% call
%   matRad_calcQualityIndicators(d,cst,refGy,refVol)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   refGy: (optional)   array of dose values used for V_XGy calculation
%                       default is [40 50 60]
%   refVol:(optional)   array of volumes (0-100) used for D_X calculation
%                       default is [2 5 95 98]
%                       NOTE: Call either both or none!
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
    refVol = [2 5 98 95];
end

numOfVois = size(cst,1);

% calculate QIs per VOI
for runVoi = 1:numOfVois
    
    indices     = cst{runVoi,4};
    numOfVoxels = numel(indices);
    
    voiPrint = sprintf('%3d %20s',cst{runVoi,1},cst{runVoi,2}); %String that will print quality indicators
    % get Dose, dose is sorted to simplify calculations
    if sum(strcmp(fieldnames(result),'RBExDose')) > 0
        relevantDose = result.RBExDose;
        doseInVoi   = sort(result.RBExDose(indices));
    else
        relevantDose = result.physicalDose;
        doseInVoi   = sort(result.physicalDose(indices));
    end
        
    % easy stats
    QI(runVoi).mean = mean(doseInVoi);
    QI(runVoi).std = std(doseInVoi);
    QI(runVoi).max = doseInVoi(end);
    QI(runVoi).min = doseInVoi(1);
    
    voiPrint = sprintf('%s - Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
                       voiPrint,QI(runVoi).mean,QI(runVoi).std,QI(runVoi).max,QI(runVoi).min,' ');
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
    
    % create VX and DX struct fieldnames at runtime and fill
    for runDX = 1:numel(refVol)
        QI(runVoi).(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));
        voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
    end
    voiPrint = sprintf('%s\n%27s',voiPrint,' ');
    for runVX = 1:numel(refGy)
        QI(runVoi).(strcat('V',num2str(refGy(runVX)))) = VX(refGy(runVX));
        voiPrint = sprintf('%sV%dGy = %6.2f%%, ',voiPrint,refGy(runVX),VX(refGy(runVX))*100);
    end
    voiPrint = sprintf('%s\n%27s',voiPrint,' ');
    
    % if current voi is a target -> calculate homogeneity and conformity
    if strcmp(cst{runVoi,3},'TARGET') > 0      
        
        % loop over target objectives and get the lowest dose objective 
        referenceDose = inf;
        for runObjective = 1:numel(cst{runVoi,6})
           % check if this is an objective that penalizes underdosing 
           if strcmp(cst{runVoi,6}(runObjective).type,'square deviation') > 0 || strcmp(cst{runVoi,6}(runObjective).type,'square underdosing') > 0
               referenceDose = min(cst{runVoi,6}(runObjective).parameter(2),referenceDose);
           end            
        end
        
        if referenceDose == inf 
            voiPrint = sprintf('%s%s',voiPrint,'Warning: target has no objective that penalizes underdosage, ');
        else
            % Conformity Index, fieldname contains reference dose
            VTarget95 = sum(doseInVoi >= 0.95*referenceDose); % number of target voxels recieving dose >= 0.95 dPres
            VTreated95 = sum(relevantDose(:) >= 0.95*referenceDose);  %number of all voxels recieving dose >= 0.95 dPres ("treated volume")
            QI(runVoi).(strcat('CI',num2str(referenceDose))) = VTarget95^2/(numOfVoxels * VTreated95); 
        
            % Homogeneity Index (one out of many), fieldname contains reference dose        
            QI(runVoi).(strcat('HI',num2str(referenceDose))) = (DX(5) - DX(95))/referenceDose * 100;
            
            voiPrint = sprintf('%sCI = %6.4f, HI = %5.2f for reference dose of %3.1f Gy\n',voiPrint,...
                               QI(runVoi).(strcat('CI',num2str(referenceDose))),QI(runVoi).(strcat('HI',num2str(referenceDose))),referenceDose);
        end 
    end
    
    fprintf('%s\n',voiPrint);    
    
end
 
result.QI = QI;

end    

