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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
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
                matRad_cfg = MatRad_Config.instance();
                handleParent = figure(...
                    'MenuBar','none',...
                    'Units','pixels',...
                    'Position',[500 500 1000 125],...
                    'Color',matRad_cfg.gui.backgroundColor,...
                    'Name','MatRad Logo',...
                    'IntegerHandle','off',...
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

            filepath = matRad_cfg.matRadSrcRoot;

            h1 = this.widgetHandle;
            
            %matRad logo
            h2 = axes(...
                'Parent',h1,...
                'Units','normalized',...
                'FontName','CMU Serif',...
                'Position',[0 0.0 0.5 0.875],...'Position',[- 0.304874274661509 -0.12225992317542 0.994397163120567 1.048719590268886],...
                'SortMethod','childorder',...
                'Tag','axesLogo');

            [im,alpha] = matRad_getLogo();
            f = image(im,'Parent',h2);
            axis(h2,'image','off');
            set(f, 'AlphaData', alpha);


            %dkfz logo
            h7 = axes(...
                'Parent',h1,...
                'Units','normalized',...
                'FontName','CMU Serif',...
                'Position',[0.5 0.0 0.5 0.85],...
                'SortMethod','childorder',...
                'Tag','axesDKFZ');

            [im,alpha] = matRad_getLogoDKFZ();

            f = image(im,'Parent',h7);
            axis(h7,'image','off');
            axis(h7,'tight');
            set(f, 'AlphaData', alpha);

            uicontrol('Parent',h1,...
                'Style','Text',...
                'Units','Normalized',...
                'Position',[0 0.85 0.5 0.15],...
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


