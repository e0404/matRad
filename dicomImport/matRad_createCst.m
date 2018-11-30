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

cst = cell(size(structures,2),6);

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
        cst{i,6}(1).type       = 'square deviation';
        cst{i,6}(1).penalty    = 800;
        cst{i,6}(1).dose       = 30;
        cst{i,6}(1).EUD        = NaN;
        cst{i,6}(1).volume     = NaN;
        cst{i,6}(1).robustness = 'none';
        
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
end
