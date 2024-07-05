function [doseCube_physicalDose, doseCube_physicalDoseRelError] = matRad_doseCubeFromDij(dij)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate dose cube from dij sparse matrix
%
% Call:
%   [doseCube.physicalDose, doseCube.physicalDoseRelError]  = matRad_doseCubeFromDij(dij)
%
% Outout:
%   doseCube.physicalDose & doseCube.physicalDose_relError
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get total physical dose
doseCube_physicalDose_dummy = zeros(prod(dij.doseGrid.dimensions),1);
doseCube_physicalDose = zeros(dij.doseGrid.dimensions);

for i = 1:size(dij.physicalDose{1,1},2)
    doseCube_physicalDose_dummy = doseCube_physicalDose_dummy + ...
        full(dij.physicalDose{1,1}(:,i));
    max(doseCube_physicalDose_dummy)
    disp(i)
end
doseCube_physicalDose(:) =  doseCube_physicalDose_dummy;

% Get error of total physical dose
doseCube_physicalDoseRelError_dummy = zeros(prod(dij.doseGrid.dimensions),1);
doseCube_physicalDoseRelError = zeros(dij.doseGrid.dimensions);

try
    for i = 1:size(dij.physicalDose{1,1},2)
        doseCube_physicalDoseRelError_dummy = doseCube_physicalDoseRelError_dummy + ...
            full(dij.physicalDose_relError{1,1}(:,i));
    end
    doseCube_physicalDoseRelError(:) =  doseCube_physicalDoseRelError_dummy;

catch
    disp('No relative error calculation available. Probably no Monte Carlo simulation run...')
end

end