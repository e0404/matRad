function matRad_applyThemeToWaitbar(hWaitbar)
% matRad_applyThemeToDlg is a helper function to apply a theme to a waitbar
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
    set(hWaitbar,'Color',matRad_cfg.gui.backgroundColor);
    txtObj = findall(hWaitbar,'type','text');
    txtObj = txtObj(1);
    if matRad_cfg.isOctave                
        patchObj = findobj(hWaitbar,'type','patch');
        set(patchObj,'FaceColor',matRad_cfg.gui.highlightColor);
        set(patchObj,'EdgeColor',matRad_cfg.gui.textColor);
        axesObj = findobj(hWaitbar,'type','axes');
        set(axesObj,'Color',matRad_cfg.gui.elementColor);
    else
        %no longer a patch object, we need java
        %jcolor = [int32(rgb * 255) 0];
        hJava = findall(hWaitbar,'type','hgjavacomponent');
        hIndic = findall(hWaitbar,'type','uiprogressindicator');
        if ~isempty(hJava)
            jcolorbg = num2cell(matRad_cfg.gui.elementColor);
            jcolorfg = num2cell(matRad_cfg.gui.highlightColor);
               
            hJava.JavaPeer.setForeground(java.awt.Color(jcolorfg{:}))
            hJava.JavaPeer.setBackground(java.awt.Color(jcolorbg{:}));
        elseif ~isempty(hIndic)
            hWaitbar.Children.Color = matRad_cfg.gui.elementColor; 
            hIndic.ProgressColor = matRad_cfg.gui.highlightColor;
        else
            %Do nothing
        end
    end 

    set(txtObj,'Color',matRad_cfg.gui.textColor,'Interpreter','none');
catch
    matRad_cfg.dispWarning('Theme could not be applied to dialog!');
end
