classdef matRad_InfoWidget_uiwrapper < matRad_Widget
% Class that defines and deploys the infoWidget to display matRad information 
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
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
        function this = matRad_InfoWidget_uiwrapper(handleParent)
            if nargin < 1                
                 handleParent = matRad_uiwrapper('uifigure',...
                     'Units','pixels',...%'characters',...
                     'Position',[150 450 190 130],...%[100 45 25 7],...
                     ...%'Visible',1,...%'on',...
                     'Color',[0.501960784313725 0.501960784313725 0.501960784313725],... 
                     ...%'IntegerHandle',0,...%'off',...
                     'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
...%                     'MenuBar','none',...
                     'Name','MatRad Info',...
                     ...%'NumberTitle',0,...%'off',... 
                     'HandleVisibility','callback',...
                     'Tag','figure1');
            end
            this = this@matRad_Widget(handleParent);
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
             h94 = this.widgetHandle;
            
            h95 = matRad_uiwrapper('uibutton',...
                 'Parent',h94,...
                 'Position',[60   2   70    30],...%[0.238095238095238 0.134831460674157 0.563492063492063 0.280898876404494],...
                 ...%'Units','pixels',...
                 'Tag','btnAbout',...
                 'FontSize',7,...
                 'text','About',...
                 'ButtonPushedFcn',@(hObject,eventdata) btnAbout_Callback(this,hObject,eventdata),...
                 'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontWeight','bold');
          
            h96 = matRad_uiwrapper('uilabel',...
                'Parent',h94,...
                ...%'Units','normalized',...
                ...%'String','v3.0.0',...
                'Position',[20   80   150    40],...%[0.227106227106227 0.752808988764045 0.523809523809524 0.191011235955056],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Tag','text15',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h97 = matRad_uiwrapper('uilabel',...
                'Parent',h94,...
                ...%'Units','normalized',...
                ...%'String','github.com/e0404/matRad',...
                'Position',[20   40   150    20],...%[0.0384615384615385 0.528089887640449 0.942307692307693 0.168539325842697],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Tag','text31',...
                'FontSize',8,...
                'FontWeight','bold' );
             
            this.createHandles();
            handles=this.handles;
            %Alter matRad Version string positioning
            vString = matRad_version();
            vPos = get(handles.text15,'Position');
            urlPos = get(handles.text31,'Position');
            btnPos = get(handles.btnAbout,'Position');
            
            %vPos([1 3]) = urlPos([1 3]);
            %vPos([1 3]) = [0 1];
%             vPos(4) = vPos(4)*1.25;
%             btnPos(2) = 0.05;
%             urlPos(2) = btnPos(2)+btnPos(4)+0.05;
%             vPos(2) = urlPos(2) + urlPos(4) + 0.05;
%             vPos(4) = 0.98 - vPos(2);
            
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
            
            [env,envver]  = matRad_getEnvironment();
            msg{end+1} = sprintf('Environment: %s v%s %s',env,envver,version('-release'));
            
            msg{end+1} = 'Web: www.matrad.org';
            msg{end+1} = 'E-Mail: contact@matrad.org';
            
            msgbox(msg,'About matRad');
            
            
            this.handles = handles;
        end
    end
end
