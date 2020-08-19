classdef matRad_PlanWidget < matRad_Widget
    
    properties
        State = false
        Machines
        Optimizations
    end
    
    properties (Constant)
        Modalities = {'photons','protons','carbon'};
    end
    
    methods
        function this = matRad_PlanWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[100 45 85 15],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],... 
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Plan',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);
            
%             if evalin('base','exist(''pln'')')
%               getPlnFromWorkspace(this);
%             else
%               setPlnDefaultValues(this);
%             end
            update(this);
            
            handles=this.handles;
            matRad_cfg = MatRad_Config.instance();
            if matRad_cfg.eduMode
                %Visisbility in Educational Mode
                eduHideHandles =   {handles.radiobutton3Dconf,...
                    handles.btnRunDAO};
                eduDisableHandles = {handles.editCouchAngle,handles.popUpMachine};
                cellfun(@(h) set(h,'Visible','Off'),eduHideHandles);
                cellfun(@(h) set(h,'Enable','Off'),eduDisableHandles);
            end
            this.handles=handles;

        end
        
        function this = initialize(this)
        end
        
        function this = update(this,evt)          
            doUpdate = true;
            if nargin == 2
                %At pln changes and at cst/cst (for Isocenter and new settings) 
                %we need to update
                doUpdate = this.checkUpdateNecessary({'pln','ct','cst'},evt);
            end
            
            if doUpdate
                if evalin('base','exist(''pln'')')
                  getPlnFromWorkspace(this);
                else
                  setPlnDefaultValues(this);
                end
            end
        end
        
        function changeWorkspace(obj)
            [env, ~] = matRad_getEnvironment();
            % handle environment
            switch env
                case 'MATLAB'
                    %the PlanWidget only changes the pln
                    evt = matRad_WorkspaceChangedEvent('pln');
                    notify(obj, 'workspaceChanged',evt);
                case 'OCTAVE'
                    matRad_notifyOctave(obj, 'workspaceChanged');
            end
        end
    end
    
    methods(Access = protected)
        function this = createLayout(this)
            h12 = this.widgetHandle;
            
            h13 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','bixel width in [mm]',...
                'Style','text',...
                'Position',[0.02 0.859315589353612 0.2 0.0950570342205324],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Tag','txtBixelWidth',...
                'UserData',[],...
                'FontSize',8,...
                'FontName','Helvetica');
            
            h14 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','5',...
                'Style','edit',...
                'Position',[0.219794344473008 0.889733840304182 0.161953727506427 0.0836501901140684],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Children',[],...
                'Tag','editBixelWidth',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h15 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Gantry Angle in °',...
                'Style','text',...
                'Position',[0.032133676092545 0.752851711026616 0.176092544987147 0.0950570342205324],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Children',[],...
                'FontSize',8,...
                'Tag','txtGantryAngle' );
            
            h16 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','0',...
                'Style','edit',...
                'Position',[0.219794344473008 0.779467680608365 0.161953727506427 0.0836501901140684],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata)standardCallback(this,hObject,eventdata),...
                'Tag','editGantryAngle',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h17 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Couch Angle in °',...
                'Style','text',...
                'Position',[0.0347043701799486 0.64638783269962 0.173521850899743 0.0950570342205324],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtCouchAngle');
            
            h18 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','0',...
                'Style','edit',...
                'Position',[0.219794344473008 0.669201520912547 0.161953727506427 0.0836501901140685],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Tag','editCouchAngle',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h19 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String',this.Modalities,...,...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.219794344473008 0.52851711026616 0.161953727506427 0.114068441064639],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata)popupRadMode_Callback(this,hObject,eventdata),...
                'Tag','popupRadMode',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h20 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Radiation Mode',...
                'Style','text',...
                'Position',[0.02 0.539923954372624 0.2 0.0950570342205324],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtRadMode');
            
            h21 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','# Fractions',...
                'Style','text',...
                'Position',[0.02 0.209125475285171 0.2 0.106463878326996],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtNumOfFractions');
            
            h22 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','30',...
                'Style','edit',...
                'Position',[0.219794344473008 0.228136882129278 0.161953727506427 0.0836501901140684],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Tag','editFraction',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h23 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','IsoCenter in [mm]',...
                'Style','text',...
                'Position',[0.02 0.330798479087452 0.2 0.091254752851711],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Tag','txtIsoCenter',...
                'FontSize',8,...
                'UserData',[],...
                'FontName','Helvetica');
            
            h24 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','0 0 0',...
                'Style','edit',...
                'Position',[0.219794344473008 0.338403041825095 0.161953727506427 0.0836501901140684],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Children',[],...
                'Enable','off',...
                'Tag','editIsoCenter',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h25 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Auto.',...
                'Style','checkbox',...
                'Value',1,...
                'Position',[0.38560411311054 0.338403041825095 0.1 0.091254752851711],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','checkIsoCenter');
            
              h26 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Run Sequencing',...
                'Style','radiobutton',...
                'Position',[0.553984575835475 0.628020880324805 0.2 0.140684410646388],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'FontSize',8,...
                'Enable','off',...
                'Tag','btnRunSequencing');
            
            h27 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Run Direct Aperture Optimization',...
                'Style','radiobutton',...
                'Position',[0.553984575835475 0.32003608945028 0.4 0.140684410646388],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata)standardCallback(this,hObject,eventdata),...
                'FontSize',8,...
                'Enable','off',...
                'Tag','btnRunDAO' );
            
            h28 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Stratification Levels',...
                'Style','text',...
                'Position',[0.553984575835475 0.502545595153702 0.2 0.102661596958175],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtSequencing' );
            
            h29 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','7',...
                'Style','edit',...
                'Position',[0.58611825192802 0.449313655990204 0.0668380462724936 0.0836501901140685],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata)standardCallback(this,hObject,eventdata),...
                'Enable','off',...
                'FontSize',8,...
                'Tag','editSequencingLevel',...
                'FontWeight','bold');
            
            h30 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String',{'Generic','generic_MCsquare'},...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.219794344473008 0.418250950570342 0.161953727506427 0.114068441064639],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) popUpMachine_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','popUpMachine',...
                'FontWeight','bold');
           
            h31 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Machine',...
                'Style','text',...
                'Position',[0.02 0.433460076045627 0.2 0.0950570342205323],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtMachine' );
            
             h32 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Set Tissue',...
                'Position',[0.401028277634961 0.110266159695817 0.109254498714653 0.0874524714828897],...
                'BackgroundColor',[0.8 0.8 0.8],...
                'Callback',@(hObject,eventdata) btnSetTissue_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'FontSize',8,...
                'Tag','btnSetTissue'); 
            
            h33 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String',{  'none'; 'const_RBExD'; 'LEMIV_effect'; 'LEMIV_RBExD' },...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.219794344473008 0.0760456273764259 0.165809768637532 0.11787072243346],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) popMenuBioOpt_Callback(this,hObject,eventdata),...
                'Tag','popMenuBioOpt',...
                'Enable', 'off',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h34 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Type of optimization',...
                'Style','text',...
                'Position',[0.0102827763496144 0.0988593155893536 0.201799485861183 0.091254752851711],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Interruptible','off',...
                'Tag','text38',...
                'FontSize',8,...
                'FontName','Helvetica' );
            
            h35 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','Dose Calculation: ',...
                'Style','text',...
                'Position',[0.5332245532245534 0.209125475285171 0.201799485861183 0.091254752851711],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Interruptible','off',...
                'Tag','text39',...
                'FontSize',8,...
                'FontName','Helvetica' );
            
            h36 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','3D conformal',...
                'Style','radiobutton',...
                'Position',[0.553224553224553 0.757869249394673 0.212121212121212 0.0847457627118644],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...                
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Enable','off',...
                'FontSize',8,...
                'Tag','radiobutton3Dconf' );
            
              h37 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','5',...
                'Style','edit',...
                'Position',[0.553224553224553 0.0760456273764259 0.09 0.11787072243346],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Children',[],...
                'Tag','editDoseX',...
                'FontSize',8,...
                'FontWeight','bold');
            
              h38 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','5',...
                'Style','edit',...
                'Position',[0.653224553224553 0.0760456273764259 0.09 0.11787072243346],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Children',[],...
                'Tag','editDoseY',...
                'FontSize',8,...
                'FontWeight','bold');
            
              h39 = uicontrol(...
                'Parent',h12,...
                'Units','normalized',...
                'String','5',...
                'Style','edit',...
                'Position',[0.753224553224553 0.0760456273764259 0.09 0.11787072243346],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) standardCallback(this,hObject,eventdata),...
                'Children',[],...
                'Tag','editDoseZ',...
                'FontSize',8,...
                'FontWeight','bold');
            
            this.createHandles();
             
        end
       
        function this = setPlnDefaultValues(this)
            
            handles = this.handles;
            
            this.getMachines()

            %
            vChar = get(handles.editGantryAngle,'String');
            if strcmp(vChar(1,1),'0') && length(vChar)==6
              set(handles.editGantryAngle,'String','0');
            end
            vChar = get(handles.editCouchAngle,'String');
            if strcmp(vChar(1,1),'0') && length(vChar)==3
              set(handles.editCouchAngle,'String','0')
            end
            
            % do not calculate / suggest isoCenter new by default
            %this.checkIsoCenter_Callback(handles.checkIsoCenter);
            set(handles.editIsoCenter,'Enable','on');

            this.handles=handles;
            updatePlnInWorkspace(this);
        end
        
        function this = getPlnFromWorkspace(this)
            pln = evalin('base', 'pln');
            handles = this.handles;
            
            % sanity check of isoCenter
            if size(pln.propStf.isoCenter,1) ~= pln.propStf.numOfBeams && size(pln.propStf.isoCenter,1) == 1
                pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * pln.propStf.isoCenter(1,:);
            elseif size(pln.propStf.isoCenter,1) ~= pln.propStf.numOfBeams && size(pln.propStf.isoCenter,1) ~= 1
                error('Isocenter in plan file are incosistent.');
            end
            
            set(handles.editBixelWidth,'String',num2str(pln.propStf.bixelWidth));
            set(handles.editGantryAngle,'String',num2str(pln.propStf.gantryAngles));
            set(handles.editCouchAngle,'String',num2str(pln.propStf.couchAngles));
            
            modIx = find(strcmp(pln.radiationMode,this.Modalities));
            set(handles.popupRadMode,'Value',modIx);
            
            getMachines(this);
            modIy = find(strcmp(pln.machine,this.Machines{modIx})); 
            set(handles.popUpMachine,'Value',modIy); 
            
            if isfield(pln.propStf,'isoCenter')
                if size(unique(pln.propStf.isoCenter,'rows'),1) == 1
                    set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.propStf.isoCenter(1,:)*10))./10), '\s+', ' '));
                    set(handles.checkIsoCenter,'Enable','on');
                    if get(handles.checkIsoCenter,'Value')
                        set(handles.editIsoCenter,'Enable','off');
                    else
                        set(handles.editIsoCenter,'Enable','on');
                    end
                    
                else
                    set(handles.editIsoCenter,'String','multiple isoCenter');
                    set(handles.editIsoCenter,'Enable','off');
                    set(handles.checkIsoCenter,'Value',0);
                    set(handles.checkIsoCenter,'Enable','off');
                end
            end
            
            set(handles.editFraction,'String',num2str(pln.numOfFractions));
            
            
            contentPopUp = get(handles.popMenuBioOpt,'String');
            ix = find(strcmp(pln.propOpt.bioOptimization,contentPopUp));
            set(handles.popMenuBioOpt,'Value',ix);
            
            set(handles.btnRunSequencing,'Value',pln.propOpt.runSequencing);
            set(handles.btnRunDAO,'Value',pln.propOpt.runDAO);
            set(handles.radiobutton3Dconf,'Value',pln.propOpt.conf3D);
            
            set(handles.editDoseX,'String',num2str(pln.propDoseCalc.doseGrid.resolution.x));
            set(handles.editDoseY,'String',num2str(pln.propDoseCalc.doseGrid.resolution.y));
            set(handles.editDoseZ,'String',num2str(pln.propDoseCalc.doseGrid.resolution.z));

            this.handles=handles;
            this.switchEnables();
        end
        
        %Update the workspace pln from the Widget
        function updatePlnInWorkspace(this)
            this.getMachines();
            handles = this.handles;

            % evalin pln (if existant) in order to decide whether isoCenter should be calculated
            % automatically
            if evalin('base','exist(''pln'',''var'')')
                pln = evalin('base','pln');
            end
            
            pln.propStf.bixelWidth      = this.parseStringAsNum(get(handles.editBixelWidth,'String'),false); % [mm] / also corresponds to lateral spot spacing for particles
            pln.propStf.gantryAngles    = this.parseStringAsNum(get(handles.editGantryAngle,'String'),true); % [???]
            pln.propStf.couchAngles     = this.parseStringAsNum(get(handles.editCouchAngle,'String'),true); % [???]
            pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
            pln.propStf.isoCenter       = this.parseStringAsNum(get(handles.editIsoCenter,'String'),true);
            
            % switch machines depending on radmode selection
            selectedMachine                     = get(handles.popUpMachine,'Value');
            popupMachines                       = get(handles.popUpMachine,'String');
            pln.machine                         = popupMachines{selectedMachine};
            

            pln.propDoseCalc.doseGrid.resolution.x = this.parseStringAsNum(get(handles.editDoseX,'String'),false);
            pln.propDoseCalc.doseGrid.resolution.y = this.parseStringAsNum(get(handles.editDoseY,'String'),false);
            pln.propDoseCalc.doseGrid.resolution.z = this.parseStringAsNum(get(handles.editDoseZ,'String'),false);
                  
            try
                ct = evalin('base','ct');
                pln.numOfVoxels     = prod(ct.cubeDim);
                pln.voxelDimensions = ct.cubeDim;
            catch
            end
            pln.numOfFractions  = this.parseStringAsNum(get(handles.editFraction,'String'),false);
            contents            = get(handles.popupRadMode,'String');
            pln.radiationMode   = contents{get(handles.popupRadMode,'Value')}; % either photons / protons / carbon
            contents            = get(handles.popUpMachine,'String');
            
            if (~strcmp(pln.radiationMode,'photons'))
                contentBioOpt = get(handles.popMenuBioOpt,'String');
                pln.propOpt.bioOptimization = contentBioOpt{get(handles.popMenuBioOpt,'Value'),:};
            else
                pln.propOpt.bioOptimization = 'none';
            end
            
            pln.propOpt.runSequencing = logical(get(handles.btnRunSequencing,'Value'));
            pln.propOpt.runDAO = logical(get(handles.btnRunDAO,'Value'));
            pln.propOpt.conf3D = logical(get(handles.radiobutton3Dconf,'Value'));
            
            
            % checkIsoCenter checkbox
            W = evalin('base','whos');
            doesPlnExist = ismember('pln',{W(:).name}) && evalin('base','exist(''cst'')') && evalin('base','exist(''ct'')');
            
            if get(handles.checkIsoCenter,'Value') && doesPlnExist
                try
                    %pln = evalin('base','pln');
                    if ~isfield(pln.propStf,'isoCenter')
                        pln.propStf.isoCenter = NaN;
                    end
                    tmpIsoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
                    if ~isequal(tmpIsoCenter,pln.propStf.isoCenter)
                        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
                        %handles.State = 1;
                        %UpdateState(handles);
                    end
                    set(handles.editIsoCenter,'String',regexprep(num2str((round(tmpIsoCenter*10))./10), '\s+', ' '));
                    set(handles.editIsoCenter,'Enable','off')
                    assignin('base','pln',pln);
                catch ME
                    warning(ME.identifier,'couldn''t set isocenter in pln update! Reason: %s\n',ME.message)
                end
            else
                set(handles.editIsoCenter,'Enable','on')
            end
            
            
            % editIsoCenter textbox
            tmpIsoCenter = str2num(get(handles.editIsoCenter,'String'));
            
            if length(tmpIsoCenter) == 3
                if sum(any(unique(pln.propStf.isoCenter,'rows')~=tmpIsoCenter))
                    pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
                    %handles.State = 1;
                    %UpdateState(handles);
                end
            else
                handles = showError(this,'EditIsoCenterCallback: Could not set iso center');
            end
            
            if evalin('base','exist(''cst'')')
                try
                    cst = evalin('base','cst');
                    if (sum(strcmp('TARGET',cst(:,3))) > 0 && get(handles.checkIsoCenter,'Value')) || ...
                            (sum(strcmp('TARGET',cst(:,3))) > 0 && ~isfield(pln.propStf,'isoCenter'))
                        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct);
                        set(handles.checkIsoCenter,'Value',1);
                    else
                        if ~strcmp(get(handles.editIsoCenter,'String'),'multiple isoCenter')
                            pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * str2num(get(handles.editIsoCenter,'String'));
                        end
                    end
                catch ME
                    warning(ME.identifier,'couldn''t set isocenter in pln update! Reason: %s\n',ME.message)
                end
            end
            
            handles.pln = pln;
            assignin('base','pln',pln);
            this.handles = handles;
            %switchEnables(this);
            changeWorkspace(this);                     
        end
    end
    
    methods(Access = private)
        function standardCallback(this, hObject, eventdata)
            %handles = this.handles;
            updatePlnInWorkspace(this);
            
            %this.handles = handles;
        end
        
        function switchEnables(this)
            handles = this.handles;
            hObject = handles.popupRadMode;
            
            contents      = cellstr(get(hObject,'String'));
            RadIdentifier = contents{get(hObject,'Value')};
            contentPopUp  = get(handles.popMenuBioOpt,'String');
            
            switch RadIdentifier
                case 'photons'
                    
                    set(handles.popMenuBioOpt,'Enable','off');
                    ix = find(strcmp(contentPopUp,'none'));
                    set(handles.popMenuBioOpt,'Value',ix);
                    set(handles.btnSetTissue,'Enable','off');
                    
                    set(handles.btnRunSequencing,'Enable','on');
                    set(handles.btnRunDAO,'Enable','on');
                    set(handles.radiobutton3Dconf,'Enable','on');
                    set(handles.txtSequencing,'Enable','on');
                    set(handles.editSequencingLevel,'Enable','on');                                           
                    
                case 'protons'                    
                    set(handles.popMenuBioOpt,'Enable','on');
                    set(handles.btnSetTissue,'Enable','on');
                    
                    set(handles.btnRunSequencing,'Enable','off');
                    set(handles.btnRunDAO,'Enable','off');
                    set(handles.radiobutton3Dconf,'Enable','off');
                    set(handles.txtSequencing,'Enable','off');
                    set(handles.editSequencingLevel,'Enable','off');
                    
                case 'carbon'
                    
                    set(handles.popMenuBioOpt,'Enable','on');
                    set(handles.btnSetTissue,'Enable','on');
                    
                    set(handles.btnRunSequencing,'Enable','off');
                    set(handles.btnRunDAO,'Enable','off');
                    set(handles.radiobutton3Dconf,'Enable','off');
                    set(handles.txtSequencing,'Enable','off');
                    set(handles.editSequencingLevel,'Enable','off');                    
            end
            
            selectedBioOpt = get(handles.popMenuBioOpt,'Value');
            contentPopUp = get(handles.popMenuBioOpt,'String');
            if strcmp(contentPopUp{selectedBioOpt},'none')
                set(handles.btnSetTissue,'Enable','off');
            else
                set(handles.btnSetTissue,'Enable','on');
            end
            this.handles = handles;
        end
        
        function manageRadModeSpecificDisplay(this)
            handles = this.handles;
            hObject = this.popupRadMode('hObject');
            
            this.handles = handles;
        end
        
        function popupRadMode_Callback(this, hObject, eventdata)
            handles = this.handles;
            contents      = cellstr(get(hObject,'String'));
            RadIdentifier = contents{get(hObject,'Value')};
            contentPopUp  = get(handles.popMenuBioOpt,'String');
            
            switch RadIdentifier
                case 'protons'
                    ix = find(strcmp(contentPopUp,'const_RBExD'));
                    set(handles.popMenuBioOpt,'Value',ix);
                    
                case 'carbon'
                    ix = find(strcmp(contentPopUp,'LEMIV_RBExD'));
                    set(handles.popMenuBioOpt,'Value',ix);
            end
            
            pln = evalin('base','pln');
            
            % new radiation modality is photons -> just keep physicalDose
            if strcmp(contents(get(hObject,'Value')),'photons')
                try
                    AllVarNames = evalin('base','who');
                    if  ismember('resultGUI',AllVarNames)
                        resultGUI = evalin('base','resultGUI');
                        if isfield(resultGUI,'alpha');    resultGUI = rmfield(resultGUI,'alpha');   end
                        if isfield(resultGUI,'beta');     resultGUI = rmfield(resultGUI,'beta');    end
                        if isfield(resultGUI,'RBExDose'); resultGUI = rmfield(resultGUI,'RBExDose');end
                        if isfield(resultGUI,'RBE');      resultGUI = rmfield(resultGUI,'RBE');     end
                        assignin('base','resultGUI',resultGUI);
                        %handles = updateIsoDoseLineCache(handles);
                    end
                catch
                end
            elseif strcmp(contents(get(hObject,'Value')),'protons')
                try
                    AllVarNames = evalin('base','who');
                    if  ismember('resultGUI',AllVarNames)
                        resultGUI = evalin('base','resultGUI');
                        if isfield(resultGUI,'alpha'); resultGUI = rmfield(resultGUI,'alpha');end
                        if isfield(resultGUI,'beta');  resultGUI = rmfield(resultGUI,'beta'); end
                        if isfield(resultGUI,'RBE');   resultGUI = rmfield(resultGUI,'RBE');  end
                        assignin('base','resultGUI',resultGUI);
                        %handles = updateIsoDoseLineCache(handles);
                    end
                catch
                end
            end
            this.handles = handles;
            updatePlnInWorkspace(this);
        end
       
%         function editIsoCenter_Callback(this, hObject, eventdata)
%             
%             handles = this.handles;
%             
%             
%             pln = evalin('base','pln');
%             tmpIsoCenter = str2num(get(hObject,'String'));
%             
%             if length(tmpIsoCenter) == 3
%                 if sum(any(unique(pln.propStf.isoCenter,'rows')~=tmpIsoCenter))
%                     pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
%                     handles.State = 1;
%                     %UpdateState(handles);
%                 end
%             else
%                 handles = showError(this,'EditIsoCenterCallback: Could not set iso center');
%             end
%             
%             assignin('base','pln',pln);
%             
%             updatePlnInWorkspace(this);
%             this.handles = handles;
%            
%         end
        
%         function checkIsoCenter_Callback(this, hObject, eventdata)
%             handles = this.handles;
%             
%             W = evalin('base','whos');
%             doesPlnExist = ismember('pln',{W(:).name});
%             
%             if get(hObject,'Value') && doesPlnExist
%                 try
%                     pln = evalin('base','pln');
%                     if ~isfield(pln.propStf,'isoCenter')
%                         pln.propStf.isoCenter = NaN;
%                     end
%                     tmpIsoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
%                     if ~isequal(tmpIsoCenter,pln.propStf.isoCenter)
%                         pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
%                         handles.State = 1;
%                         %UpdateState(handles);
%                     end
%                     set(handles.editIsoCenter,'String',regexprep(num2str((round(tmpIsoCenter*10))./10), '\s+', ' '));
%                     set(handles.editIsoCenter,'Enable','off')
%                     assignin('base','pln',pln);
%                 catch ME
%                     warning(ME.identifier,'couldn''t set isocenter in pln update! Reason: %s\n',ME.message)
%                 end
%             else
%                 set(handles.editIsoCenter,'Enable','on')
%             end
%             
%             updatePlnInWorkspace(this); 
%             this.handles = handles;
%         end
        
        function popUpMachine_Callback(this, hObject, eventdata)
            % M�GLICHER FEHLER WEGEN VALUE WERT!
            handles = this.handles;
             contents = cellstr(get(hObject,'String'));
             MachineIdentifier = contents{get(hObject,'Value')};
            % contentPopUp = get(handles.)
            flag=checkRadiationComposition(this);
            if ~flag
                this.showWarning(['No base data available for machine: ' MachineIdentifier '. Selecting default machine.']);
                set(handles.popUpMachine,'Value',1);
            end
            getMachines(this);
            pln = evalin('base','pln');
            
            % M�GLICHEE FEHLER HIER VALUE UND GENERIC WERDEN VERGLICHEN
            if strcmp(contents(get(hObject,'Value')),'Generic')
                try
                    AllVarNames = evalin('base','who');
                    if  ismember('resultGUI',AllVarNames)
                        resultGUI = evalin('base','resultGUI');
                        if isfield(resultGUI,'alpha');    resultGUI = rmfield(resultGUI,'alpha');   end
                        if isfield(resultGUI,'beta');     resultGUI = rmfield(resultGUI,'beta');    end
                        if isfield(resultGUI,'RBExDose'); resultGUI = rmfield(resultGUI,'RBExDose');end
                        if isfield(resultGUI,'RBE');      resultGUI = rmfield(resultGUI,'RBE');     end
                        assignin('base','resultGUI',resultGUI);
                        %handles = updateIsoDoseLineCache(handles);
                    end
                catch
                end
            % M�GLICHEE FEHLER HIER VALUE UND GENERIC WERDEN VERGLICHEN
            elseif strcmp(contents(get(hObject,'Value')),'generic_MCsquare')
                try
                    AllVarNames = evalin('base','who');
                    if  ismember('resultGUI',AllVarNames)
                        resultGUI = evalin('base','resultGUI');
                        if isfield(resultGUI,'alpha'); resultGUI = rmfield(resultGUI,'alpha');end
                        if isfield(resultGUI,'beta');  resultGUI = rmfield(resultGUI,'beta'); end
                        if isfield(resultGUI,'RBE');   resultGUI = rmfield(resultGUI,'RBE');  end
                        assignin('base','resultGUI',resultGUI);
                        %handles = updateIsoDoseLineCache(handles);
                    end
                catch
                end
            end
               
            this.handles = handles;
            updatePlnInWorkspace(this); 
        end
        
        function btnSetTissue_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            if evalin('base','exist(''cst'')') && evalin('base','exist(''pln'')') 
                try
                    %parse variables from base-workspace
                    cst = evalin('base','cst');
                    pln = evalin('base','pln');
                    
                    
                    fileName = [pln.radiationMode '_' pln.machine];
                    load(fileName);
                    
                    % check for available cell types characterized by alphaX and betaX
                    for i = 1:size(machine.data(1).alphaX,2)
                        CellType{i} = [num2str(machine.data(1).alphaX(i)) ' ' num2str(machine.data(1).betaX(i))];
                    end
                    
                    %fill table data array
                    for i = 1:size(cst,1)
                        data{i,1} = cst{i,2};
                        data{i,2} = [num2str(cst{i,5}.alphaX) ' ' num2str(cst{i,5}.betaX)];
                        data{i,3} = (cst{i,5}.alphaX / cst{i,5}.betaX );
                    end
                    
                    Width  = 400;
                    Height = 200 + 20*size(data,1);
                    ScreenSize = get(0,'ScreenSize');
                    % show "set tissue parameter" window
                    figHandles = get(0,'Children');
                    if ~isempty(figHandles)
                        IdxHandle = strcmp(get(figHandles,'Name'),'Set Tissue Parameters');
                    else
                        IdxHandle = [];
                    end
                    
                    %check if window is already exists
                    if any(IdxHandle)
                        IdxTable = find(strcmp({figHandles(IdxHandle).Children.Type},'uitable'));
                        set(figHandles(IdxHandle).Children(IdxTable), 'Data', []);
                        figTissue = figHandles(IdxHandle);
                        %set focus
                        figure(figTissue);
                    else
                        figTissue = figure('Name','Set Tissue Parameters','Color',[.5 .5 .5],'NumberTitle','off','Position',...
                            [ceil(ScreenSize(3)/2) ceil(ScreenSize(4)/2) Width Height]);
                    end
                    
                    % define the tissue parameter table
                    cNames = {'VOI','alphaX betaX','alpha beta ratio'};
                    columnformat = {'char',CellType,'numeric'};
                    
                    tissueTable = uitable('Parent', figTissue,'Data', data,'ColumnEditable',[false true false],...
                        'ColumnName',cNames, 'ColumnFormat',columnformat,'Position',[50 150 10 10]);
                    set(tissueTable,'CellEditCallback',@(hObject,eventdata) tissueTable_CellEditCallback(this,hObject,eventdata));
                    % set width and height
                    currTablePos = get(tissueTable,'Position');
                    currTableExt = get(tissueTable,'Extent');
                    currTablePos(3) = currTableExt(3);
                    currTablePos(4) = currTableExt(4);
                    set(tissueTable,'Position',currTablePos);
                    
                    % define two buttons with callbacks
                    uicontrol('Parent', figTissue,'Style', 'pushbutton', 'String', 'Save&Close',...
                        'Position', [Width-(0.25*Width) 0.1 * Height 70 30],...
                        'Callback', @(hpb,eventdata)SaveTissueParameters(this,hpb,eventdata));
                    
                    uicontrol('Parent', figTissue,'Style', 'pushbutton', 'String', 'Cancel&Close',...
                        'Position', [Width-(0.5*Width) 0.1 * Height 80 30],...
                        'Callback', 'close');
                catch ME
                    warning(ME.identifier,'couldn''t set isocenter in pln update! Reason: %s\n',ME.message)
                end
            end
            this.handles = handles;
            %updatePlnInWorkspace(this);
        end
        
        function popMenuBioOpt_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            pln = evalin('base','pln');
            contentBioOpt = get(handles.popMenuBioOpt,'String');
            NewBioOptimization = contentBioOpt(get(handles.popMenuBioOpt,'Value'),:);
            
                if (strcmp(pln.propOpt.bioOptimization,'LEMIV_effect') && strcmp(NewBioOptimization,'LEMIV_RBExD')) ||...
                        (strcmp(pln.propOpt.bioOptimization,'LEMIV_RBExD') && strcmp(NewBioOptimization,'LEMIV_effect'))
                    % do nothing - re-optimization is still possible
                elseif ((strcmp(pln.propOpt.bioOptimization,'const_RBE') && strcmp(NewBioOptimization,'none')) ||...
                        (strcmp(pln.propOpt.bioOptimization,'none') && strcmp(NewBioOptimization,'const_RBE'))) && isequal(pln.radiationMode,'protons')
                    % do nothing - re-optimization is still possible
                end
                
            this.handles = handles;
            updatePlnInWorkspace(this);
        end
                   
        function getMachines(this)
            matRad_cfg = MatRad_Config.instance();
            %seach for availabes machines
            handles = this.handles;
            this.Machines=cell(size(this.Modalities));
            %Loop over all modalities to find machine per modalitiy
            for i = 1:length(this.Modalities)
                pattern = [this.Modalities{1,i} '_*'];
                if isdeployed
                    baseroot = [ctfroot filesep 'matRad'];
                else
                    baseroot = matRad_cfg.matRadRoot;
                end
                Files = dir([baseroot filesep 'basedata' filesep pattern]);
                
                for j = 1:length(Files)
                    if ~isempty(Files)
                        MachineName = Files(j).name(numel(this.Modalities{1,i})+2:end-4);
                        this.Machines{i}{j} = MachineName;
%                         if isfield(handles,'Machines')
%                             if sum(strcmp(handles.Machines,MachineName)) == 0
%                                 handles.Machines{size(handles.Machines,2)+1} = MachineName;
%                             end
%                         else
%                             handles.Machines = cell(1);
%                             handles.Machines{1} = MachineName;
%                         end
                    end
                end
            end
            
            selectedRadMod = get(handles.popupRadMode,'Value');
            nMachines = numel(this.Machines{selectedRadMod});
            selectedMachine = get(handles.popUpMachine,'Value');
            
            if get(handles.popUpMachine,'Value') > nMachines
                selectedMachine = 1;
            end
                            
            set(handles.popUpMachine,'Value',selectedMachine,'String',this.Machines{selectedRadMod});
            this.handles = handles;
        end
        
        function number = parseStringAsNum(this,stringIn,isVector)
            if isnumeric(stringIn)
                number = stringIn;
            else
                number = str2num(stringIn);
                if isempty(number) || length(number) > 1 && ~isVector
                    warndlg(['could not parse all parameters (pln, optimization parameter)']);
                    number = NaN;
                elseif isVector && iscolumn(number)
                    number = number';
                end
            end
        end
        
        function flag = checkRadiationComposition(this)
            handles = this.handles;
            
            flag = true;
            contents = cellstr(get(handles.popUpMachine,'String'));
            Machine = contents{get(handles.popUpMachine,'Value')};
            contents = cellstr(get(handles.popupRadMode,'String'));
            radMod = contents{get(handles.popupRadMode,'Value')};
            
            if isdeployed
                baseroot = [ctfroot filesep 'matRad'];
            else
                baseroot = [fileparts(mfilename('fullpath')) filesep '..'];
            end
            FoundFile = dir([baseroot filesep 'basedata' filesep  radMod '_' Machine '.mat']);
            
%             if isdeployed
%                 FoundFile = dir([ctfroot filesep 'matRad' filesep radMod '_' Machine '.mat']);
%             else
%                 FoundFile = dir([fileparts(mfilename('fullpath')) filesep  '..' filesep radMod '_' Machine '.mat']);
%             end
            if isempty(FoundFile)
              %  this.showWarning(['No base data available for machine: ' Machine '. Selecting default machine.']);
                flag = false;
              %  set(handles.popUpMachine,'Value',1);
            end
            this.handles = handles;
        end
        
        function SaveTissueParameters(this,~, ~)
            cst = evalin('base','cst');
            % get handle to uiTable
            figHandles = get(0,'Children');
            IdxHandle  = find(strcmp(get(figHandles,'Name'),'Set Tissue Parameters'));
            % find table in window
            
            figHandleChildren = get(figHandles(IdxHandle),'Children');
            IdxTable   = find(strcmp(get(figHandleChildren,'Type'),'uitable'));
            uiTable    = figHandleChildren(IdxTable);
            % retrieve data from uitable
            data       = get(uiTable,'data');
            
            for i = 1:size(cst,1)
                for j = 1:size(data,1)
                    if strcmp(cst{i,2},data{j,1})
                        alphaXBetaX = str2num(data{j,2});
                        cst{i,5}.alphaX = alphaXBetaX(1);
                        cst{i,5}.betaX  = alphaXBetaX(2);
                    end
                end
            end
            assignin('base','cst',cst);
            close
            %handles.State = 2;
            %UpdateState(handles);
            
            %this.handles = handles;
            updatePlnInWorkspace(this); 
        end
        
        function tissueTable_CellEditCallback(this,hObject, eventdata)
            if eventdata.Indices(2) == 2
                alphaXBetaX = str2num(eventdata.NewData);
                data = get(hObject,'Data');
                data{eventdata.Indices(1),3} = alphaXBetaX(1)/alphaXBetaX(2);
                set(hObject,'Data',data);
            end
        end
    end
end
