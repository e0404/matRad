function [outputArg1,outputArg2] = matRad_createMCsquareBaseDataFile(filename,machine,focusIx)
% matRad - create base data from requested proton machine for MCsquare
% 
% call
%   matRad_createMCsquareBaseDataFile(machine)
%
% input
%   filename:   name of the output file
%   machine:    matRad machine (protons)
%   focusIx:    Focus setting (MCsquare only allows one focus setting)
%
%
% output
%
% References
%   -
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~strcmp(machine.meta.radiationMode,'protons')
    %matRad_dispToConsole(
    error('MCsquare only working with proton base data');    
end

%nozzleToIso = machine.meta.BAMStoIsoDist;

nEnergies = numel(machine.data);

SAD = machine.meta.SAD;

fileID = fopen(filename,'w');

%Header
%fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
fprintf(fileID,'# %s\n',machine.meta.description);
fprintf(fileID,'# created by %s on %s\n\n',machine.meta.created_by,machine.meta.created_on);

nozzleToIso = 500;
fprintf(fileID,'Nozzle exit to Isocenter distance\n');
fprintf(fileID,'%.1f\n\n',nozzleToIso);

fprintf(fileID,'SMX to Isocenter distance\n');
fprintf(fileID,'%.1f\n\n',SAD);

fprintf(fileID,'SMY to Isocenter distance\n');
fprintf(fileID,'%.1f\n\n',SAD);

fprintf(fileID,'Beam parameters\n%d energies\n\n',nEnergies);

NominalEnergy = [machine.data(:).energy]';
rowSz = size(NominalEnergy);
MeanEnergy = NominalEnergy; %Needs a correct value
EnergySpread = 0.01*NominalEnergy; %Just an assumption of 1% energy accuracy

ProtonsMU = ones(rowSz)*1e6; %Needs a correct value

%Calculate width & divergence at iso in sigma
fwhmIso = machine.meta.LUT_bxWidthminFWHM(2,focusIx);
sigmaIso = 0.5* fwhmIso / sqrt(2*log(2));
spotDiv = sigmaIso / SAD;

fclose(fileID);

if ~isfield(machine.data,'sigma')
    Weight1 = ones(rowSz);%1 - [machine.data(:).weight]';
    SpotSize1x = ones(rowSz)*sigmaIso;
    Divergence1x = ones(rowSz)*spotDiv;
    Correlation1x = ones(rowSz)*0;
    SpotSize1y = ones(rowSz)*sigmaIso;
    Divergence1y = ones(rowSz)*spotDiv;
    Correlation1y = ones(rowSz)*0;
    
    Weight2 = zeros(rowSz);%[machine.data(:).weight]';
    SpotSize2x = ones(rowSz)*sigmaIso;
    Divergence2x = ones(rowSz)*spotDiv;
    Correlation2x = ones(rowSz)*0;
    SpotSize2y = ones(rowSz)*sigmaIso;
    Divergence2y = ones(rowSz)*spotDiv;
    Correlation2y = ones(rowSz)*0;
else
    Weight1 = ones(rowSz);
    SpotSize1x = ones(rowSz)*sigmaIso;
    Divergence1x = ones(rowSz)*spotDiv;
    Correlation1x = ones(rowSz)*0;
    SpotSize1y = ones(rowSz)*sigmaIso;
    Divergence1y = ones(rowSz)*spotDiv;
    Correlation1y = ones(rowSz)*0;
    
    Weight2 = ones(rowSz)*0;
    SpotSize2x = ones(rowSz)*sigmaIso;
    Divergence2x = ones(rowSz)*spotDiv;
    Correlation2x = ones(rowSz)*0;
    SpotSize2y = ones(rowSz)*sigmaIso;
    Divergence2y = ones(rowSz)*spotDiv;
    Correlation2y = ones(rowSz)*0;
end

dataTable = table(NominalEnergy,MeanEnergy,EnergySpread,ProtonsMU,Weight1,SpotSize1x,Divergence1x,Correlation1x,SpotSize1y,Divergence1y,Correlation1y,Weight2,SpotSize2x,Divergence2x,Correlation2x,SpotSize2y,Divergence2y,Correlation2y);

writetable(dataTable,['tmp_' filename],'Delimiter','\t');

fileID = fopen(['tmp_' filename],'r');
tableTxt = fread(fileID);
fclose(fileID);

fileID = fopen(filename,'a');
fwrite(fileID,tableTxt);
fclose(fileID);

delete(['tmp_' filename]);

end


