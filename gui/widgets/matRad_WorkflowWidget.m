classdef matRad_WorkflowWidget < matRad_Widget
    
    properties
    end
    
    methods
        function this = matRad_WorkflowWidget(handleParent)
            if nargin < 1
                matRad_cfg = MatRad_Config.instance();
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[130 45 150 15],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Workflow',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
                set(handleParent,'Units','normalized')
                pos = get(handleParent,'Position');
                pos(1:2) = [0.5 0.5] - pos(3:4)./2;
                set(handleParent,'Position');
                
            end
            this = this@matRad_Widget(handleParent);
        end
        
        function this = initialize(this)
            this.update();
        end
        
        function this = update(this,evt)
            getFromWorkspace(this);
            %updateInWorkspace(this);
        end       
        
        
        % moved so it can be called from the toolbar button
        % H74 Callback
        function btnLoadMat_Callback(this, hObject, event)
            handles = this.handles;
            [FileName, FilePath] = uigetfile('*.mat');
            if FileName == 0 % user pressed cancel --> do nothing.
                return;
            end

            try
                % delete existing workspace - parse variables from base workspace
                AllVarNames = evalin('base','who');
                RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
                
                for i = 1:length(RefVarNames)
                    if sum(ismember(AllVarNames,RefVarNames{i}))>0
                        evalin('base',['clear ', RefVarNames{i}]);
                    end
                end
                
                % read new data
                load([FilePath FileName]);
                
            catch ME
                this.handles=handles;
                getFromWorkspace(this);
                showError(this,'LoadMatFileFnc: Could not load *.mat file',ME);                
                return
            end
            
            try
                %cst = generateCstTable(this,cst);
                %handles.TableChanged = false;
                %set(handles.popupTypeOfPlot,'Value',1);
                %cst = matRad_computeVoiContoursWrapper(cst,ct);
                
                assignin('base','ct',ct);
                assignin('base','cst',cst);
                
            catch ME
                showError(this,'LoadMatFileFnc: Could not load *.mat file',ME);
            end
            
            % check if a optimized plan was loaded
            if exist('stf','var')
                assignin('base','stf',stf);
            end
            if exist('pln','var')
                assignin('base','pln',pln);
            end
            if exist('dij','var')
                assignin('base','dij',dij);
            end
            
            if exist('resultGUI','var')
                assignin('base','resultGUI',resultGUI);
            end
            
            this.handles=handles;
            %updateInWorkspace(this);
            this.changedWorkspace();
            %getFromWorkspace(this); %update the buttons
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
                       
            parent = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
                   
            h72 = this.addControlToGrid([2 5],...
                'Style','text',...
                'String','Status:',...                                
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Tag','txtStatus',...
                'FontSize',round(matRad_cfg.gui.fontSize*1.2));
            
            h73 = this.addControlToGrid([3 5],...
                'String','no data loaded',...
                'Style','text',...               
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Tag','txtInfo',...
                'FontSize',round(matRad_cfg.gui.fontSize*1.2));
            
            hMatLoad = this.addControlToGrid([2 1],...
                'String','Load  *.mat data',...                
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnLoadMat_Callback(this,hObject,eventdata),...
                'Tag','btnLoadMat');
            
             hDijCalc = this.addControlToGrid([3 1],...
                'String','Calc. Dose Influence',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnCalcDose_Callback(this,hObject,eventdata),...
                'Tag','btnCalcDose');
            
            hOpt = this.addControlToGrid([4 1],...
                'Parent',parent,...                
                'String','Optimize',...                
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnOptimize_Callback(this,hObject,eventdata),...
                'Tag','btnOptimize');
            
            hLoadDicom = this.addControlToGrid([2 2],...
                'String','Load DICOM',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnLoadDicom_Callback(this,hObject,eventdata),...
                'Tag','btnLoadDicom');
      
            
            hRefresh = this.addControlToGrid([1 1],...
                'String','Refresh',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnRefresh_Callback(this,hObject,eventdata),...
                'Tag','btnRefresh');
            
            hRecalc = this.addControlToGrid([4 2],...
                'String','Recalculate Dose',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) pushbutton_recalc_Callback(this,hObject,eventdata),...
                'Tag','pushbutton_recalc');
            
            hKeep = this.addControlToGrid([5 1],...
                'String','Save/Keep Result',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btnSaveToGUI_Callback(this,hObject,eventdata),...
                'Tag','btnSaveToGUI');
            
            hExportBin = this.addControlToGrid([5 2],...
                'String','Export Binary',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) btn_export_Callback(this,hObject,eventdata),...
                'Children',[],...
                'Tag','btn_export');
            
            hImportDose = this.addControlToGrid([4 3],...
                'String','Import Dose',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) importDoseButton_Callback(this, hObject,eventdata),...
                'Tag','importDoseButton');
            
            hImportBin = this.addControlToGrid([2 3],...
                'String','Import from Binary',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) pushbutton_importFromBinary_Callback(this,hObject,eventdata),...
                'TooltipString','Imports a patient data set from binary datafiles describing CT and segmentations',...
                'Tag','pushbutton_importFromBinary');
            
            hExportDicom = this.addControlToGrid([5 3],...
                'String','Export Dicom',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'Callback',@(hObject,eventdata) exportDicomButton_Callback(this, hObject,eventdata),...
                'Tag','exportDicomButton');
            
            this.createHandles();
            
             handles=this.handles;
            matRad_cfg = MatRad_Config.instance();
            if matRad_cfg.eduMode
                %Visisbility in Educational Mode
                eduHideHandles =   {handles.pushbutton_importFromBinary,...
                    handles.btnLoadDicom,...
                    handles.btn_export,...
                    handles.exportDicomButton,...
                    handles.importDoseButton};
                cellfun(@(h) set(h,'Visible','Off'),eduHideHandles);
            end
            this.handles=handles;
        end
        
        function this = getFromWorkspace(this)
            handles = this.handles;
            
            % no data loaded, disable the buttons
            set(handles.txtInfo,'String','no data loaded');
            set(handles.btnCalcDose,'Enable','off');
            set(handles.btnOptimize ,'Enable','off');
            set(handles.pushbutton_recalc,'Enable','off');
            set(handles.btnSaveToGUI,'Enable','off');
            set(handles.importDoseButton,'Enable','off');
            set(handles.btn_export,'Enable','off');
            set(handles.exportDicomButton,'Enable','off');
            
            if evalin('base','exist(''pln'')')
                
                if evalin('base','exist(''ct'')') && ...
                        evalin('base','exist(''cst'')')
                    
                    % ct cst and pln available; ready for dose calculation
                    set(handles.txtInfo,'String','ready for dose calculation');
                    set(handles.btnCalcDose,'Enable','on');
                    set(handles.btn_export,'Enable','on');
                    set(handles.exportDicomButton,'Enable','on');
                    
                    if evalin('base','exist(''resultGUI'')')
                        % plan is optimized
                        % check if dij, stf and pln match                        
                        if matRad_comparePlnDijStf(evalin('base','pln'),evalin('base','stf'),evalin('base','dij'))
                            set(handles.txtInfo,'String','plan is optimized');
                            set(handles.btnOptimize ,'Enable','on'); 
                        end
                        
                        set(handles.pushbutton_recalc,'Enable','on');
                        set(handles.btnSaveToGUI,'Enable','on');
                        % resultGUI struct needs to be available to import dose
                        % otherwise inconsistent states can be achieved
                        set(handles.importDoseButton,'Enable','on');
                        
                    elseif evalin('base','exist(''dij'')') &&  evalin('base','exist(''stf'')')
                        % check if dij, stf and pln match
                        if matRad_comparePlnDijStf(evalin('base','pln'),evalin('base','stf'),evalin('base','dij'))
                            % plan is ready for optimization
                            set(handles.txtInfo,'String','ready for optimization');
                            set(handles.btnOptimize ,'Enable','on');
                        end
                    end
                end
            end
            this.handles=handles;
        end
        
    end
    methods (Access = private)
        
        function h = addControlToGrid(this,gridPos,varargin)
            matRad_cfg = MatRad_Config.instance();
            parent = this.widgetHandle;
            
            %Use a 5 x 5 grid
            pos = this.computeGridPos(gridPos,[5 5]);                  
            
            h = uicontrol('Parent',parent,...
                'Units','normalized',...
                'Position',pos,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                varargin{:});
        end
        
        
        % H75 Callback
        function btnCalcDose_Callback(this, hObject, eventdata)
            % hObject    handle to btnCalcDose (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % http://stackoverflow.com/questions/24703962/trigger-celleditcallback-before-button-callback
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/332613
            % wait some time until the CallEditCallback is finished
            % Callback triggers the cst saving mechanism from GUI
            
            handles = this.handles;
            try
                % indicate that matRad is busy
                % change mouse pointer to hour glass
                Figures = gcf;%findobj('type','figure');
                set(Figures, 'pointer', 'watch');
                drawnow;
                % disable all active objects
                InterfaceObj = findobj(Figures,'Enable','on');
                set(InterfaceObj,'Enable','off');
                
                
                % read plan from gui and save it to workspace
                %handles=getPlnFromGUI(this);
                
                % get default iso center as center of gravity of all targets if not
                % already defined
                pln = evalin('base','pln');
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'CalcDoseCallback: Error in preprocessing!',ME);
                return;
            end
            
            % generate steering file
            try
                currPln = evalin('base','pln');
                %                     % if we run 3d conf opt -> hijack runDao to trigger computation of
                %                     % connected bixels
                %                     if strcmp(pln.radiationMode,'photons') && get(handles.radiobutton3Dconf,'Value')
                %                         currpln.propOpt.runDAO = true;
                %                     end
                stf = matRad_generateStf(evalin('base','ct'),...
                    evalin('base','cst'),...
                    currPln);
                assignin('base','stf',stf);
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'CalcDoseCallback: Error in steering file generation!',ME);
                return;
            end
            
            % carry out dose calculation
            try
                if strcmp(pln.radiationMode,'photons')
                    dij = matRad_calcPhotonDose(evalin('base','ct'),stf,pln,evalin('base','cst'));
                elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
                    dij = matRad_calcParticleDose(evalin('base','ct'),stf,pln,evalin('base','cst'));
                end
                
                % assign results to base worksapce
                assignin('base','dij',dij);
                
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'CalcDoseCallback: Error in dose calculation!',ME);
                return;
            end
            
            % change state from busy to normal
            set(Figures, 'pointer', 'arrow');
            set(InterfaceObj,'Enable','on');
            this.handles = handles;
            this.changedWorkspace('stf','dij');
            %getFromWorkspace(this);
        end
        
        % H76 Callback
        function btnOptimize_Callback(this, hObject, eventdata)
            % hObject    handle to btnOptimize (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles = this.handles;
            try
                % indicate that matRad is busy
                % change mouse pointer to hour glass
                Figures = gcf;%findobj('type','figure');
                set(Figures, 'pointer', 'watch');
                drawnow;
                % disable all active objects
                InterfaceObj = findobj(Figures,'Enable','on');
                set(InterfaceObj,'Enable','off');
                
                pln = evalin('base','pln');
                ct  = evalin('base','ct');
                
                % optimize
                [resultGUIcurrentRun,usedOptimizer] = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),pln);
                if pln.propOpt.conf3D && strcmp(pln.radiationMode,'photons')
                    resultGUIcurrentRun.w = resultGUIcurrentRun.w * ones(dij.totalNumOfBixels,1);
                    resultGUIcurrentRun.wUnsequenced = resultGUIcurrentRun.w;
                end
                
                %if resultGUI already exists then overwrite the "standard" fields
                AllVarNames = evalin('base','who');
                if  ismember('resultGUI',AllVarNames)
                    resultGUI = evalin('base','resultGUI');
                    sNames = fieldnames(resultGUIcurrentRun);
                    oldNames = fieldnames(resultGUI);
                    if(length(oldNames) > length(sNames))
                        for j = 1:length(oldNames)
                            if strfind(oldNames{j}, 'beam')
                                resultGUI = rmfield(resultGUI, oldNames{j});
                            end
                        end
                    end
                    for j = 1:length(sNames)
                        resultGUI.(sNames{j}) = resultGUIcurrentRun.(sNames{j});
                    end
                else
                    resultGUI = resultGUIcurrentRun;
                end
                assignin('base','resultGUI',resultGUI);
                %
                %                     % set some values
                %                     if handles.plane == 1
                %                         set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,2)/ct.resolution.x));
                %                     elseif handles.plane == 2
                %                         set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,1)/ct.resolution.y));
                %                     elseif handles.plane == 3
                %                         set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,3)/ct.resolution.z));
                %                     end
                
                %                     handles.State = 3;
                %                     handles.SelectedDisplayOptionIdx = 1;
                %                     if strcmp(pln.radiationMode,'carbon') || (strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD'))
                %                         handles.SelectedDisplayOption = 'RBExDose';
                %                     else
                %                         handles.SelectedDisplayOption = 'physicalDose';
                %                     end
                %                     handles.selectedBeam = 1;
                % check IPOPT status and return message for GUI user if no DAO or
                % particles
                if ~pln.propOpt.runDAO || ~strcmp(pln.radiationMode,'photons')
                    CheckOptimizerStatus(this,usedOptimizer,'Fluence')
                end
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'OptimizeCallback: Could not optimize!',ME);
                return;
            end
            
            % perform sequencing and DAO
            try
                
                %% sequencing
                if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
                    resultGUI = matRad_siochiLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
                        ,str2double(get(handles.editSequencingLevel,'String')));
                    
                    assignin('base','resultGUI',resultGUI);
                end
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'OptimizeCallback: Could not perform sequencing',ME);
                return;
            end
            
            try
                %% DAO
                if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
                    handles = showWarning(handles,['Observe: You are running direct aperture optimization' filesep 'This is experimental code that has not been thoroughly debugged - especially in combination with constrained optimization.']);
                    [resultGUI,usedOptimizer] = matRad_directApertureOptimization(evalin('base','dij'),evalin('base','cst'),...
                        resultGUI.apertureInfo,resultGUI,pln);
                    assignin('base','resultGUI',resultGUI);
                    % check IPOPT status and return message for GUI user
                    CheckOptimizerStatus(this,usedOptimizer,'DAO');
                end
                
                if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
                    matRad_visApertureInfo(resultGUI.apertureInfo);
                end
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'OptimizeCallback: Could not perform direct aperture optimization',ME);
                return;
            end
            
            % change state from busy to normal
            set(Figures, 'pointer', 'arrow');
            set(InterfaceObj,'Enable','on');
            %                 handles.dispWindow{3,1}  = [];   % reset dose ranges
            %                 handles.dispWindow{3,2}  = [];   % reset min max dose values
            %                 handles.rememberCurrAxes = false;
            %                 handles.IsoDose.Levels   = 0;  % ensure to use default iso dose line spacing
            %                 handles.cBarChanged      = true;
            %                 handles = updateIsoDoseLineCache(handles);
            %                 UpdatePlot(handles);
            %                 handles.rememberCurrAxes = true;
            this.handles = handles;
            
            this.changedWorkspace('resultGUI');
            %getFromWorkspace(this);
        end
        
        % H77 Callback
        function btnLoadDicom_Callback(this, hObject, event)
            % hObject    handle to btnLoadDicom (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles = this.handles;
            try
                % delete existing workspace - parse variables from base workspace
                %                     set(handles.popupDisplayOption,'String','no option available');
                AllVarNames = evalin('base','who');
                RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
                for i = 1:length(RefVarNames)
                    if sum(ismember(AllVarNames,RefVarNames{i}))>0
                        evalin('base',['clear ', RefVarNames{i}]);
                    end
                end
                matRad_importDicomWidget;
                
            catch ME
                showError(this,'DicomImport: Could not import data', ME);
            end
            
            this.handles = handles;
            %changeWorkspace(this);
            %getFromWorkspace(this);
        end
        
        % H78 Callback - button: refresh
        function btnRefresh_Callback(this, hObject, event)
            % notify so all widgets refresh
            this.changedWorkspace();
            %getFromWorkspace(this);
            
            %                 handles = this.handles;
            %
            %                 %% parse variables from base workspace
            %                 AllVarNames = evalin('base','who');
            %                 try
            %                     if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
            %                         ct  = evalin('base','ct');
            %                         cst = evalin('base','cst');
            %                         cst=generateCstTable(this,cst);
            %                         cst = matRad_computeVoiContoursWrapper(cst,ct);
            %                         assignin('base','cst',cst);
            %
            %                         changeWorkspace(this);
            %
            %                     elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
            %                         handles = showError(this,'GUI OpeningFunc: could not find cst file');
            %                     elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
            %                         handles = showError(this,'GUI OpeningFunc: could not find ct file');
            %                     end
            %                 catch ME
            %                    handles = showError(this,'GUI OpeningFunc: Could not load ct and cst file. Reason: ', ME);
            %                 end
            %
            %                 this.handles = handles;
        end
        
        
        % H79 Callback
        function pushbutton_recalc_Callback(this, hObject, eventdata)
            
            handles = this.handles;
            
            try
                % indicate that matRad is busy
                % change mouse pointer to hour glass
                Figures = gcf;%findobj('type','figure');
                set(Figures, 'pointer', 'watch');
                drawnow;
                % disable all active objects
                InterfaceObj = findobj(Figures,'Enable','on');
                set(InterfaceObj,'Enable','off');
                
                % get all data from workspace
                pln       = evalin('base','pln');
                stf       = evalin('base','stf');
                ct        = evalin('base','ct');
                cst       = evalin('base','cst');
                resultGUI = evalin('base','resultGUI');
                
                %                         % get weights of the selected cube
                %                         Content = get(handles.popupDisplayOption,'String');
                %                         SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
                %                         Suffix = strsplit(SelectedCube,'_');
                %                         if length(Suffix)>1
                %                             Suffix = ['_' Suffix{2}];
                %                         else
                %                             Suffix = '';
                %                         end
                
                if sum([stf.totalNumOfBixels]) ~= length(resultGUI.w)%(['w' Suffix]))
                    warndlg('weight vector does not corresponding to current steering file');
                    return
                end
                
                % change isocenter if that was changed and do _not_ recreate steering
                % information
                for i = 1:numel(pln.propStf.gantryAngles)
                    stf(i).isoCenter = pln.propStf.isoCenter(i,:);
                end
                
                % recalculate influence matrix
                if strcmp(pln.radiationMode,'photons')
                    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
                elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
                    dij = matRad_calcParticleDose(ct,stf,pln,cst);
                end
                
                % recalculate cubes in resultGUI
                resultGUIreCalc = matRad_calcCubes(resultGUI.w,dij,cst); %(['w' Suffix])
                
                % delete old variables to avoid confusion
                if isfield(resultGUI,'effect')
                    resultGUI = rmfield(resultGUI,'effect');
                    resultGUI = rmfield(resultGUI,'RBExDose');
                    resultGUI = rmfield(resultGUI,'RBE');
                    resultGUI = rmfield(resultGUI,'alpha');
                    resultGUI = rmfield(resultGUI,'beta');
                end
                
                % overwrite the "standard" fields
                sNames = fieldnames(resultGUIreCalc);
                for j = 1:length(sNames)
                    resultGUI.(sNames{j}) = resultGUIreCalc.(sNames{j});
                end
                
                % assign results to base worksapce
                assignin('base','dij',dij);
                assignin('base','resultGUI',resultGUI);
                
                
                %                         % show physicalDose of newly computed state
                %                         handles.SelectedDisplayOption = 'physicalDose';
                %                         set(handles.popupDisplayOption,'Value',find(strcmp('physicalDose',Content)));
                
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                
                %                         handles.cBarChanged = true;
                %                         handles = updateIsoDoseLineCache(handles);
                
                this.handles = handles;
                this.workspaceChanged('dij','resultGUI');
                
            catch ME
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                showError(this,'CalcDoseCallback: Error in dose recalculation!',ME);
                return;
                
            end
            
        end
        
        % H80 Callback
        function btnSaveToGUI_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            Width  = 400;
            Height = 200;
            ScreenSize = get(0,'ScreenSize');
            
            % show "Provide result name" window
            figHandles = get(0,'Children');
            if ~isempty(figHandles)
                IdxHandle = strcmp(get(figHandles,'Name'),'Provide result name');
            else
                IdxHandle = [];
            end
            
            %check if window is already exists
            if any(IdxHandle)
                figDialog = figHandles(IdxHandle);
                %set focus
                figure(figDialog);
            else
                figDialog = dialog('Position',[ceil(ScreenSize(3)/2) ceil(ScreenSize(4)/2) Width Height],'Name','Provide result name','Color',[0.5 0.5 0.5]);
                
                uicontrol('Parent',figDialog,...
                    'Style','text',...
                    'Position',[20 Height - (0.35*Height) 350 60],...
                    'String','Please provide a decriptive name for your optimization result:','FontSize',10,'BackgroundColor',[0.5 0.5 0.5]);
                
                uicontrol('Parent',figDialog,...
                    'Style','edit',...
                    'Position',[30 60 350 60],...
                    'String','Please enter name here...','FontSize',10,'BackgroundColor',[0.55 0.55 0.55]);
                
                uicontrol('Parent', figDialog,'Style', 'pushbutton', 'String', 'Save','FontSize',10,...
                    'Position', [0.42*Width 0.1 * Height 70 30],...
                    'Callback', @(hpb,eventdata)SaveResultToGUI(this,hpb,eventdata));
            end
            
            uiwait(figDialog);
            this.handles = handles;
            %changeWorkspace(this);
            %getFromWorkspace(this);
        end
        
        function SaveResultToGUI(this, ~, ~)
            AllFigHandles = get(0,'Children');
            ixHandle      = strcmp(get(AllFigHandles,'Name'),'Provide result name');
            uiEdit        = get(AllFigHandles(ixHandle),'Children');
            
            if strcmp(get(uiEdit(2),'String'),'Please enter name here...')
                
                formatOut = 'mmddyyHHMM';
                Suffix = ['_' datestr(now,formatOut)];
            else
                % delete special characters
                Suffix = get(uiEdit(2),'String');
                logIx = isstrprop(Suffix,'alphanum');
                Suffix = ['_' Suffix(logIx)];
            end
            
            pln       = evalin('base','pln');
            resultGUI = evalin('base','resultGUI');
            
            if isfield(resultGUI,'physicalDose')
                resultGUI.(['physicalDose' Suffix])  = resultGUI.physicalDose;
            end
            if isfield(resultGUI,'w')
                resultGUI.(['w' Suffix])             = resultGUI.w;
            end
            
            
            if ~strcmp(pln.propOpt.bioOptimization,'none')
                
                if isfield(resultGUI,'RBExDose')
                    resultGUI.(['RBExDose' Suffix]) = resultGUI.RBExDose;
                end
                
                if strcmp(pln.radiationMode,'carbon') == 1
                    if isfield(resultGUI,'effect')
                        resultGUI.(['effect' Suffix])= resultGUI.effect;
                    end
                    
                    if isfield(resultGUI,'RBE')
                        resultGUI.(['RBE' Suffix]) = resultGUI.RBE;
                    end
                    if isfield(resultGUI,'alpha')
                        resultGUI.(['alpha' Suffix]) = resultGUI.alpha;
                    end
                    if isfield(resultGUI,'beta')
                        resultGUI.(['beta' Suffix]) = resultGUI.beta;
                    end
                end
            end
            
            close(AllFigHandles(ixHandle));
            assignin('base','resultGUI',resultGUI);
            
            this.changedWorkspace('resultGUI');
            %getFromWorkspace(this);
        end
        
        % H81 Callback
        function btn_export_Callback(this, hObject, eventdata)
            % hObject    handle to btn_export (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            try
                matRad_exportWidget;
            catch ME
                showError(this,'Could not export data.  Reason: ', ME);
            end
        end
        
        % H82 Callback
        function importDoseButton_Callback(this, hObject, eventdata)
            % hObject    handle to importDoseButton (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles = this.handles;
            
            extensions{1} = '*.nrrd';
            [filenames,filepath,~] = uigetfile(extensions,'MultiSelect','on');
            
            try
                if ~iscell(filenames)
                    tmp = filenames;
                    filenames = cell(1);
                    filenames{1} = tmp;
                end
                
                ct = evalin('base','ct');
                resultGUI = evalin('base','resultGUI');
                
                for filename = filenames
                    [~,name,~] = fileparts(filename{1});
                    [cube,~] = matRad_readCube(fullfile(filepath,filename{1}));
                    if ~isequal(ct.cubeDim, size(cube))
                        errordlg('Dimensions of the imported cube do not match with ct','Import failed!','modal');
                        continue;
                    end
                    
                    fieldname = ['import_' matlab.lang.makeValidName(name, 'ReplacementStyle','delete')];
                    resultGUI.(fieldname) = cube;
                end
                
                assignin('base','resultGUI',resultGUI);
            catch ME
                this.handles = handles;
                showError(this,'Dose Import: Could not import data.  Reason: ', ME);
                return;
            end
            this.handles = handles;
            this.changedWorkspace('resultGUI');
            %getFromWorkspace(this);
        end
        
        % H83 Callback
        function pushbutton_importFromBinary_Callback(this, hObject, eventdata)
            % hObject    handle to pushbutton_importFromBinary (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles = this.handles;
            
            try
                % delete existing workspace - parse variables from base workspace
                %set(handles.popupDisplayOption,'String','no option available');
                AllVarNames = evalin('base','who');
                RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
                ChangedVarNames = {};
                for i = 1:length(RefVarNames)
                    if sum(ismember(AllVarNames,RefVarNames{i}))>0
                        evalin('base',['clear ', RefVarNames{i}]);
                    end
                end
                
                %call the gui
                h=matRad_importWidget;
                uiwait(h.widgetHandle);
                
%                 %Check if we have the variables in the workspace
%                 if evalin('base','exist(''cst'',''var'')') == 1 && evalin('base','exist(''ct'',''var'')') == 1
%                     cst = evalin('base','cst');
%                     ct = evalin('base','ct');
%                     %cst = generateCstTable(this,cst);
%                     %                         handles.TableChanged = false;
%                     %                         set(handles.popupTypeOfPlot,'Value',1);
%                     
%                     % compute HU values
%                     if ~isfield(ct, 'cubeHU')
%                         ct = matRad_electronDensitiesToHU(ct);
%                     end
%                     %                         if ~isfield(ct, 'cubeHU')
%                     %                             handles.cubeHUavailable = false;
%                     %                         else
%                     %                             handles.cubeHUavailable = true;
%                     %                         end
%                     
%                     %                         % precompute contours
%                     %                         cst = precomputeContours(this,ct,cst);
%                     
%                     assignin('base','ct',ct);
%                     assignin('base','cst',cst);
%                     
%                     if evalin('base','exist(''pln'',''var'')')
%                         assignin('base','pln',pln);
%                         %                             setPln(handles);
%                         %                         else
%                         %                             getPlnFromGUI(handles);
%                         %                             setPln(handles);
%                     end
%                     
%                    
%                     
%                 end
                this.handles = handles;
                this.changedWorkspace();
            catch ME                
                this.handles = handles;
                getFromWorkspace(this);
                showError(this,'Binary Patient Import: Could not import data.  Reason: ', ME);
                return;
            end
            
            
            %getFromWorkspace(this);
        end
        
        % H84 Callback
        function exportDicomButton_Callback(this, hObject, eventdata)
            % hObject    handle to exportDicom (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            try                
                matRad_exportDicomWidget;
            catch ME
                showError(this,'DicomImport: Could not export data', ME);
            end
        end
        
        function CheckOptimizerStatus(this, usedOptimizer,OptCase)
            
            [statusmsg,statusflag] = usedOptimizer.GetStatus();
            
            if statusflag == 0 || statusflag == 1
                status = 'none';
            else
                status = 'warn';
            end
            
            msgbox(['Optimizer finished with status ' num2str(statusflag) ' (' statusmsg ')'],'Optimizer',status,'modal');
        end
    end
end


