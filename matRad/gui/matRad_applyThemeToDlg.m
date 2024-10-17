function matRad_applyThemeToDlg(hDlgBox)
% matRad_applyThemeToDlg is a helper function to apply a theme tho dialog
% boxes
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
try
    if matRad_cfg.isOctave
        txtObj = findall(hDlgBox,'type','text');
        txtObj = txtObj(1);
        okBtn = findall(hDlgBox,'type','uicontrol','style','pushbutton');
        set(findall(hDlgBox,'type','uipanel'),'BackgroundColor',matRad_cfg.gui.backgroundColor);
    else
        txtObj = findall(hDlgBox,'tag','MessageBox');
        okBtn = findall(hDlgBox,'tag','OKButton');
        hDlgBox.Color = matRad_cfg.gui.backgroundColor;
    end
    
    set(txtObj,'Color',matRad_cfg.gui.textColor);
    set(okBtn,'BackgroundColor',matRad_cfg.gui.elementColor,'ForegroundColor',matRad_cfg.gui.textColor);
catch
    matRad_cfg.dispWarning('Theme could not be applied to dialog!');
end