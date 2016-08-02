function resultGUI = matRad_postprocessing(resultGUI, dij,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad postprosseing function accounting for
%       minimum number of particles per spot
%       (scan path)
%       ...
% 
% call
%   resultGUI = matRad_fluenceOptimization(resultGUI, dij,cst,pln)
%
% input
%   resultGUI   struct containing optimized fluence vector
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% issue warning if biological optimization impossible
if sum(strcmp(pln.bioOptimization,{'effect','RBExD'}))>0 && (~isfield(dij,'mAlphaDose') || ~isfield(dij,'mSqrtBetaDose'))
    warndlg('Alpha and beta matrices for effect based and RBE optimization not available - physical optimization is carried out instead.');
    pln.bioOptimization = 'none';
end

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

end

%%manipulate fluence vector
w = resultGUI.w;
Imin = pln.minNrParticles/1e6;
lw = length(w);
for i = 1:lw
    if(w(i) < Imin/2)
        w(i) = 0;
    elseif(w(i) > Imin/2 && w(i) < Imin)
        w(i) = Imin;        
    end
end

%%calc dose
d = matRad_backProjection(w,dij,'none');

resultGUI.finalDose = reshape(d{1},dij.dimensions);
resultGUI.finalw = w;

%%calc difference to optimized dose (not necessary, can be deleted)
relIntDoseDif = (1-sum(resultGUI.physicalDose(:))/sum(resultGUI.finalDose(:)))*100;

fprintf(['Relative difference in integral dose: ' num2str(relIntDoseDif) '%%\n']);

