function hGUI = matRadGUI(varargin)
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


persistent hMatRadGUI;


p = inputParser;
addParameter(p,'devMode',false,@(x) validateModeValue(x));
addParameter(p,'eduMode',false,@(x) validateModeValue(x));
p.KeepUnmatched = true; %No error with incorrect parameters

parse(p,varargin{:});
parsedInput = p.Results;

if ischar(parsedInput.devMode) || isstring(parsedInput.devMode)
    parsedInput.devMode = str2double(parsedInput.devMode);
end

if ischar(parsedInput.eduMode) || isstring(parsedInput.eduMode)
    parsedInput.eduMode = str2double(parsedInput.eduMode);
end

matRad_cfg = MatRad_Config.instance();

%If devMode is true, error dialogs will include the full stack trace of the error
%If false, only the basic error message is shown (works for errors that
%handle the MException object)
matRad_cfg.devMode = logical(parsedInput.devMode);
if matRad_cfg.devMode
    disp('matRadGUI starting in developer mode!');
end

%Enables simple educational mode which removes certain functionality from
%the GUI
matRad_cfg.eduMode = logical(parsedInput.eduMode);
if matRad_cfg.eduMode
    disp('matRadGUI starting in educational mode!');
end

if matRad_cfg.disableGUI
    matRad_cfg.dispInfo('matRad GUI disabled in matRad_cfg!\n');
    return;
end

handleValid = true;

try
    handleValid = ishandle(hMatRadGUI.guiHandle);
catch

    handleValid = false;
end

if handleValid
    figure(hMatRadGUI.guiHandle);
else
    matRad_rc(false);
    hMatRadGUI = matRad_MainGUI;
end

if nargout > 0
    hGUI = hMatRadGUI;
end


end


%Validates the attributes for the command line Modes
function validateModeValue(x)
%When passed from OS terminal (or inline in Command Window) everything is a string
if isdeployed || ischar(x) || isstring(x)
    x=str2double(x);

end
validateattributes(x,{'logical','numeric'},{'scalar','binary'});
end
