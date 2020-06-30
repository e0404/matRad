classdef matRad_LogoWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_LogoWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'PaperUnits','inches',...
                    'MenuBar','none',...
                    'Units','characters',...
                    'Position',[135.8 70.7692307692308 500 18.3846153846154],...
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
            
            if isdeployed
                filepath = [ctfroot filesep 'matRad'];
            else
                filepath = fileparts(mfilename('fullpath'));
            end
            
            h1 = this.widgetHandle;
            
            h2 = axes(...
                'Parent',h1,...
                'FontName','CMU Serif',...
                'Position',[-0.304874274661509 -0.12225992317542 0.994397163120567 1.048719590268886],...'Position',[- 0.304874274661509 -0.12225992317542 0.994397163120567 1.048719590268886],...
                'SortMethod','childorder',...
                'Tag','axesLogo');
            
            [im, ~, alpha] = imread([filepath filesep '..' filesep 'gfx' filesep 'matrad_logo.png']);
            f = image(im,'Parent',h2);
            axis(h2,'image','off');
            set(f, 'AlphaData', alpha);
            
                      
            h7 = axes(...
                'Parent',h1,...
                'FontName','CMU Serif',...
                'Position',[0.29482269503546 0.053969270166453 0.799187620889749 0.6017029449423816],...
                'SortMethod','childorder',...
                'Tag','axesDKFZ');
           
            [im, ~, alpha] = imread([filepath filesep '..' filesep 'gfx' filesep 'DKFZ_Logo.png']);
            f = image(im,'Parent',h7);
            axis(h7,'image','off');
            set(f, 'AlphaData', alpha);
            
            
            this.createHandles();
            
        end
    end
    
    % NO CALLBACKS
end


