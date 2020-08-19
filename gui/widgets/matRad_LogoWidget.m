classdef matRad_LogoWidget < matRad_Widget
    
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
                'Position',[0 0 0.4 1],...'Position',[- 0.304874274661509 -0.12225992317542 0.994397163120567 1.048719590268886],...
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
                'Position',[0.4 0 0.6 1],...
                'SortMethod','childorder',...
                'Tag','axesDKFZ');
           
            
            f = image(im,'Parent',h7);
            axis(h7,'image','off');
            axis(h7,'tight');
            set(f, 'AlphaData', alpha);
            
            
            this.createHandles();
            
        end
    end
    
    % NO CALLBACKS
end


