function cst = matRad_createCst(structures)
% matRad function to create a cst struct upon dicom import
% 
% call
%   cst = matRad_createCst(structures)
%
% input
%   structures:     matlab struct containing information about rt structure
%                   set (generated with matRad_importDicomRtss and 
%                   matRad_convRtssContours2Indices)
%
% output
%   cst:            matRad cst struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

nStructures = size(structures,2);
cst = cell(nStructures,6);

%Create set of default colors
defaultColors = colorcube(nStructures);

for i = 1:size(structures,2)
    cst{i,1} = i - 1; % first organ has number 0    
    cst{i,2} = structures(i).structName;
    
    if ~isempty(regexpi(cst{i,2},'tv')) || ...
       ~isempty(regexpi(cst{i,2},'target')) || ...
       ~isempty(regexpi(cst{i,2},'gtv')) || ...
       ~isempty(regexpi(cst{i,2},'ctv')) || ...
       ~isempty(regexpi(cst{i,2},'ptv')) || ...
       ~isempty(regexpi(cst{i,2},'boost')) || ...
       ~isempty(regexpi(cst{i,2},'tumor'))
        
        cst{i,3} = 'TARGET';
        
        cst{i,5}.Priority = 1;
     
        % default objectives for targets
        objective = DoseObjectives.matRad_SquaredDeviation;
        objective.penalty = 800;
        objective.parameters = {30};  %Default reference Dose
        cst{i,6}{1} = struct(objective);
        
    else
        
        cst{i,3} = 'OAR';
        
        cst{i,5}.Priority = 2;
        
        cst{i,6} = []; % define no OAR dummy objcetives   
    
    end
    
    cst{i,4}{1} = structures(i).indices;
    
    % set default parameter for biological planning
    cst{i,5}.alphaX  = 0.1;
    cst{i,5}.betaX   = 0.05;
    cst{i,5}.Visible = 1;
    if isfield(structures(i),'structColor') && ~isempty(structures(i).structColor)
        cst{i,5}.visibleColor = structures(i).structColor' ./ 255;
    else
        cst{i,5}.visibleColor = defaultColors(i,:);
        matRad_cfg.dispInfo('No color information for structure %d "%s". Assigned default color [%f %f %f]\n',i,cst{i,2},defaultColors(i,1),defaultColors(i,2),defaultColors(i,3));
    end
end
