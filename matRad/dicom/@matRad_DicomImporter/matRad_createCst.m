function obj = matRad_createCst(obj)
% matRad function to create a cst struct upon dicom import
% 
% In your object, there must be a property that contains matlab structure 
% containing information about rt structure set (generated with 
% matRad_importDicomRtss and matRad_convRtssContours2Indices)
%
% Output - matRad cst structure
%
% call
%   obj = matRad_createCst(obj)
%
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

nStructures = size(obj.importRtss.structures,2);
obj.cst = cell(nStructures,6);

%Create set of default colors
defaultColors = colorcube(nStructures);

for i = 1:size(obj.importRtss.structures,2)
    obj.cst{i,1} = i - 1; % first organ has number 0    
    obj.cst{i,2} = obj.importRtss.structures(i).structName;
    
    if ~isempty(regexpi(obj.cst{i,2},'tv', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'target', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'gtv', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'ctv', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'ptv', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'boost', 'once')) || ...
       ~isempty(regexpi(obj.cst{i,2},'tumor', 'once'))
        
        obj.cst{i,3} = 'TARGET';
        
        obj.cst{i,5}.Priority = 1;
     
        % default objectives for targets
        objective = DoseObjectives.matRad_SquaredDeviation;
        objective.penalty = 800;
        objective.parameters = {30};  %Default reference Dose
        obj.cst{i,6}{1} = struct(objective);
        
    else
        
        obj.cst{i,3} = 'OAR';
        
        obj.cst{i,5}.Priority = 2;
        
        obj.cst{i,6} = []; % define no OAR dummy objcetives   
    
    end
    
    obj.cst{i,4}{1} = obj.importRtss.structures(i).indices;
    
    % set default parameter for biological planning
    obj.cst{i,5}.alphaX  = 0.1;
    obj.cst{i,5}.betaX   = 0.05;
    obj.cst{i,5}.Visible = 1;
    if isfield(obj.importRtss.structures(i),'structColor') && ~isempty(obj.importRtss.structures(i).structColor)
        obj.cst{i,5}.visibleColor = obj.importRtss.structures(i).structColor' ./ 255;
    else
        obj.cst{i,5}.visibleColor = defaultColors(i,:);
        matRad_cfg.dispInfo('No color information for structure %d "%s". Assigned default color [%f %f %f]\n',i,obj.cst{i,2},defaultColors(i,1),defaultColors(i,2),defaultColors(i,3));
    end
end
