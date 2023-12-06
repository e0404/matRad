classdef matRad_InfoWidget < matRad_Widget
    % matRad_InfoWidget class to generate GUI widget to display system and
    % version information 
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
        function this = matRad_InfoWidget(handleParent)
            if nargin < 1
                matRad_cfg = MatRad_Config.instance();
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.45 0.45 0.1 0.1],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,... 
                    'IntegerHandle','off',...                    
                    'MenuBar','none',...
                    'Name','matRad Info',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
            end
            this = this@matRad_Widget(handleParent);
        end
    end
    
        methods (Access = protected)
        function this = createLayout(this)
            h94 = this.widgetHandle;
            matRad_cfg = MatRad_Config.instance();
            txt = sprintf('Info about\nsoftware environment & version\nmatRad version & branch');
            %About button
            h95 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                'String','About',...
                'Tooltip', txt,...
                'Position',[0.2 0.14 0.6 0.28],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnAbout_Callback(this,hObject,eventdata),...
                'Tag','btnAbout',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'FontName',matRad_cfg.gui.fontName);
            
            %Position String
            h96 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                'Style','text',...
                'Position',[0.1 0.75 0.8 0.2],...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','text15',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'FontName',matRad_cfg.gui.fontName);
            
            %URL Position String
            h97 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                'Style','text',...
                'Position',[0.05 0.5 0.9 0.17],...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','text31',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight','bold');
            
            this.createHandles();
            handles=this.handles;
            %Alter matRad Version string positioning
            vString = matRad_version();
            vPos = get(handles.text15,'Position');
            urlPos = get(handles.text31,'Position');
            btnPos = get(handles.btnAbout,'Position');
            
            %vPos([1 3]) = urlPos([1 3]);
            vPos([1 3]) = [0 1];
            vPos(4) = vPos(4)*1.25;
            btnPos(2) = 0.05;
            urlPos(2) = btnPos(2)+btnPos(4)+0.05;
            vPos(2) = urlPos(2) + urlPos(4) + 0.05;
            vPos(4) = 0.98 - vPos(2);
            
            set(handles.btnAbout,'Position',btnPos);
            set(handles.text31,'String','www.matRad.org','Position',urlPos,'Enable','inactive','ButtonDownFcn', @(~,~) web('www.matrad.org','-browser'));
            set(handles.text15,'String',vString,'Position',vPos);
            this.handles=handles;
        end
    end
    
    methods (Access = protected) 
        function btnAbout_Callback(this, hObject, event)
            handles = this.handles;
            %msgbox({'https://github.com/e0404/matRad/' 'email: matrad@dkfz.de'},'About');
            
            matRad_cfg = MatRad_Config.instance();
            [~,matRadVer] = matRad_version;
            
            msg{1} = ['matRad ''' matRadVer.name '''']; %Name
            if matRad_cfg.eduMode
                msg{1} = [msg{1} ' Educational'];
            end
            msg{end+1} = sprintf('v%d.%d.%d',matRadVer.major,matRadVer.minor,matRadVer.patch); %Version Number
            if isdeployed
                msg{end+1} = 'Standalone Version';
            elseif ~isempty(matRadVer.branch) && ~isempty(matRadVer.commitID)
                msg{end+1} = sprintf('Git: Branch %s, commit %s',matRadVer.branch,matRadVer.commitID(1:8));
            end
                        
            msg{end+1} = sprintf('Environment: %s v%s %s',matRad_cfg.env,matRad_cfg.envVersion,version('-release'));
            
            msg{end+1} = 'Web: www.matrad.org';
            msg{end+1} = 'E-Mail: contact@matrad.org';
            
            msg{end+1} = 'MATRAD IS NOT A MEDICAL PRODUCT AND THEREFORE NOT SUITABLE FOR CLINICAL USE!';
            
            msgbox(msg,'About matRad');
            
            
            this.handles = handles;
        end
    end
end
