function hGUI = matRadParetoGUI(varargin)
% matRad compatability function to call the matRad_MainGUI
%   The function checks input parameters and handles the GUI as a
%   singleton, so following calls will not create new windows
%
% call
%   hGUI = matRadGUI
%   matRadGUI
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


persistent hMatRadParetoGUI;

%Initialize navigation structures -> Check for retStruct etc.
retStruct = evalin('base','retStruct');
finds = retStruct.finds;
weights = retStruct.weights;

%Matrix storing objective function values of all calculated plans
assignin('base','finds',finds);

%Matrix storing the plans' associated weights
assignin('base','weights',weights);

matRad_ParetoGUI();







