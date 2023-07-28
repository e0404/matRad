classdef matRad_LogoWidget < matRad_Widget
    % matRad_InfoWidget class to display GUI logo widget 
    % 
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
    end
    
    methods
        function this = matRad_LogoWidget(handleParent)            
            
            if nargin < 1
                handleParent = figure(...
                    'MenuBar','none',...
                    'Units','pixels',...
                    'Position',[500 500 1000 200],...
                    'Color',[0.5 0.5 0.5],...
                    'Name','MatRad Logo',...
                    'HandleVisibility','callback',...
                    'Tag','figure_importDialog',...
                    'WindowStyle','normal',... 'PaperSize',[8.5 11],...
                    'PaperType','usletter');
            end
            this = this@matRad_Widget(handleParent);
        end
        
        
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            %mfile = which(mfilename);           
            %[filepath] = fileparts(mfile);
            
            matRad_cfg = MatRad_Config.instance();
            
            if isdeployed
                filepath = [ctfroot filesep 'matRad'];
            else
                filepath = matRad_cfg.matRadRoot;
            end
            
            h1 = this.widgetHandle;
            
            
            %matRad Logo
            [im, ~, alpha] = imread([filepath filesep 'gfx' filesep 'matrad_logo.png']);
            dim = size(im);           
            
            h2 = axes(...
                'Parent',h1,...
                'Units','normalized',...
                'FontName','CMU Serif',...
                'Position',[0 0.2 0.4 0.8],...'Position',[- 0.304874274661509 -0.12225992317542 0.994397163120567 1.048719590268886],...
                'SortMethod','childorder',...
                'Tag','axesLogo');
            
            f = image(im,'Parent',h2);
            axis(h2,'image','off');
            set(f, 'AlphaData', alpha);
            
            
            %dkfz logo
            [im, ~, alpha] = imread([filepath filesep 'gfx' filesep 'DKFZ_Logo.png']);

            h7 = axes(...
                'Parent',h1,...
                'Units','normalized',...
                'FontName','CMU Serif',...
                'Position',[0.4 0.15 0.6 0.85],...
                'SortMethod','childorder',...
                'Tag','axesDKFZ');
           
            
            f = image(im,'Parent',h7);
            axis(h7,'image','off');
            axis(h7,'tight');
            set(f, 'AlphaData', alpha);
            
            hDisc = uicontrol('Parent',h1,...
                'Style','Text',...
                'Units','Normalized',...
                'Position',[0 0.05 1 0.15],...
                'String','matRad is not a medical product! DO NOT USE IT CLINICALLY!',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...                
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            this.createHandles();
            
        end
    end
    
    % NO CALLBACKS
end


