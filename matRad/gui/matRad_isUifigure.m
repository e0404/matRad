function tf = matRad_isUifigure(hFig)
% Checks if the figure was created with uifigure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
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
            
if matRad_cfg.isOctave
    tf = false;
else  
    if verLessThan('Matlab','9.0')      %version < 2016a (release of uifigs)
        tf = false;
    elseif verLessThan('Matlab','9.5')  % 16a <= version < 2018b
        tf = ~isempty(matlab.ui.internal.dialog.DialogHelper.getFigureID(hFig));
    else                                % version >= 2018b
        tf = matlab.ui.internal.isUIFigure(hFig);
    end
end

end

