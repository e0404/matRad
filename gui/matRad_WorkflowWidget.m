classdef matRad_WorkflowWidget < matRad_Widget
    
    properties
    end
    
    methods
        function this = matRad_WorkflowWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[130 45 80 8],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],... 
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Workflow',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);
        end
        
        function this = initialize(this)
            this.update();
        end
        
        function this = update(this)
            getFromWorkspace(this);
            %updateInWorkspace(this);
        end
       
        function changeWorkspace(obj)
            notify(obj, 'workspaceChanged');
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
                    handles = showError(this,'LoadMatFileFnc: Could not load *.mat file',ME);
                    
                    this.handles=handles;   
                    getFromWorkspace(this);  
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
                    handles = showError(this,'LoadMatFileFnc: Could not load *.mat file',ME);
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
               changeWorkspace(this);
               %getFromWorkspace(this); %update the buttons
            end
    end
    
        methods (Access = protected)
            function this = createLayout(this)
                h71 = this.widgetHandle;
                
                  h72 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Status:',...
                    'Style','text',...
                    'Position',[0.318250377073907 0.107438016528926 0.120663650075415 0.12396694214876],...
                    'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                    'Tag','text13',...
                    'FontWeight','bold' );
                
                h73 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','no data loaded',...
                    'Style','text',...
                    'Position',[0.414781297134238 0.0247933884297521 0.371040723981901 0.214876033057851],...
                    'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                    'Tag','txtInfo',...
                    'FontWeight','bold' );

                  h74 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Load  *.mat data',...
                    'Position',[0.151866151866152 0.810126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnLoadMat_Callback(this,hObject,eventdata),...
                    'Tag','btnLoadMat',...
                    'FontWeight','bold'); 
                
                 h75 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Calc. influence Mx',...
                    'Position',[0.35006435006435 0.810126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnCalcDose_Callback(this,hObject,eventdata),...
                    'Tag','btnCalcDose',...
                    'FontWeight','bold');
                
                h76 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Optimize',...
                    'Position',[0.544401544401541 0.810126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnOptimize_Callback(this,hObject,eventdata),...
                    'Tag','btnOptimize',...
                    'FontWeight','bold');
                
                h77 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Load DICOM',...
                    'Position',[0.151866151866152 0.60126582278481 0.177606177606178 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnLoadDicom_Callback(this,hObject,eventdata),...
                    'Tag','btnLoadDicom',...
                    'FontWeight','bold' );
                
                h78 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Refresh',...
                    'Position',[0.0154440154440154 0.810126582278481 0.0849420849420849 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnRefresh_Callback(this,hObject,eventdata),...
                    'Tag','btnRefresh',...
                    'FontWeight','bold' );
                
                h79 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Recalc',...
                    'Position',[0.543114543114543 0.60126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) pushbutton_recalc_Callback(this,hObject,eventdata),...
                    'Tag','pushbutton_recalc',...
                    'FontWeight','bold' );
                
                h80 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Save to GUI',...
                    'Position',[0.738738738738737 0.810126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btnSaveToGUI_Callback(this,hObject,eventdata),...
                    'Tag','btnSaveToGUI',...
                    'FontWeight','bold');
                
                h81 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Export',...
                    'Position',[0.74002574002574 0.60126582278481 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) btn_export_Callback(this,hObject,eventdata),...
                    'Children',[],...
                    'Tag','btn_export',...
                    'FontWeight','bold');
                
                h82 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Import Dose',...
                    'Position',[0.738738738738738 0.392405063291139 0.178893178893179 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) importDoseButton_Callback(this, hObject,eventdata),...
                    'Tag','importDoseButton',...
                    'FontWeight','bold' );
                
                h83 = uicontrol(...
                    'Parent',h71,...
                    'Units','normalized',...
                    'String','Import from Binary',...
                    'Position',[0.151866151866152 0.392405063291139 0.177606177606178 0.145569620253165],...
                    'BackgroundColor',[0.8 0.8 0.8],...
                    'Callback',@(hObject,eventdata) pushbutton_importFromBinary_Callback(this,hObject,eventdata),...
                    'TooltipString','Imports a patient data set from binary datafiles describing CT and segmentations',...
                    'Tag','pushbutton_importFromBinary',...
                    'FontWeight','bold');
                
                this.createHandles();
                
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
                   
                if evalin('base','exist(''pln'')')
                    
                    if evalin('base','exist(''ct'')') && ...
                        evalin('base','exist(''cst'')')
                    
                        % ct cst and pln available; ready for dose calculation
                        set(handles.txtInfo,'String','ready for dose calculation');
                        set(handles.btnCalcDose,'Enable','on');
                        set(handles.btn_export,'Enable','on');
                        
                        if evalin('base','exist(''resultGUI'')')
                            
                            % plan is optimized
                            set(handles.txtInfo,'String','plan is optimized');
                            set(handles.btnOptimize ,'Enable','on');
                            set(handles.pushbutton_recalc,'Enable','on');
                            set(handles.btnSaveToGUI,'Enable','on');
                            % resultGUI struct needs to be available to import dose
                            % otherwise inconsistent states can be achieved
                            set(handles.importDoseButton,'Enable','on');
                            
                        elseif evalin('base','exist(''dij'')')
                            % plan is ready for optimization
                            set(handles.txtInfo,'String','ready for optimization');
                            set(handles.btnOptimize ,'Enable','on');
                            set(handles.btnOptimize ,'Enable','on');
                        end
                    end
                end
                this.handles=handles;
            end
            
        end
        methods (Access = private)
            
           
            
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
                    handles = showError(this,'CalcDoseCallback: Error in preprocessing!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
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
                    handles = showError(this,'CalcDoseCallback: Error in steering file generation!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
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
                    handles = showError(this,'CalcDoseCallback: Error in dose calculation!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
                    return;
                end
                
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                this.handles = handles;
                changeWorkspace(this);
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
                    % wait until the table is updated
                    %btnTableSave_Callback([],[],handles); %We don't need it?
                    
%                     % if a critical change to the cst has been made which affects the dij matrix
%                     if handles.DijCalcWarning == true
%                         
%                         choice = questdlg('Overlap priorites of OAR constraints have been edited, a new OAR VOI was added or a critical row constraint was deleted. A new Dij calculation might be necessary.', ...
%                             'Title','Cancel','Calculate Dij then Optimize','Optimze directly','Optimze directly');
%                         
%                         switch choice
%                             case 'Cancel'
%                                 set(Figures, 'pointer', 'arrow');
%                                 set(InterfaceObj,'Enable','on');
%                                 % guidata(hObject,handles);
%                                 this.handles = handles;
%                                 
%                                 return;
%                             case 'Calculate dij again and optimize'
%                                 handles.DijCalcWarning = false;
%                                 btnCalcDose_Callback(hObject, eventdata, handles)
%                             case 'Optimze directly'
%                                 handles.DijCalcWarning = false;
%                         end
%                     end
                    
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
                    handles = showError(this,'OptimizeCallback: Could not optimize!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
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
                    handles = showError(this,'OptimizeCallback: Could not perform sequencing',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
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
                    handles = showError(this,'OptimizeCallback: Could not perform direct aperture optimization',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    this.handles = handles;
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
                
               changeWorkspace(this);
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
                    handles = showError(this,'DicomImport: Could not import data', ME);
                end
                
                this.handles = handles;
                changeWorkspace(this);
                %getFromWorkspace(this);
            end
            
            % H78 Callback - button: refresh
            function btnRefresh_Callback(this, hObject, event)
                % notify so all widgets refresh
                changeWorkspace(this);
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
                        changeWorkspace(this);
                        %getFromWorkspace(this);
                        
                    catch ME
                        handles = showError(this,'CalcDoseCallback: Error in dose recalculation!',ME);
                        
                        % change state from busy to normal
                        set(Figures, 'pointer', 'arrow');
                        set(InterfaceObj,'Enable','on');
                        this.handles = handles;
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
                
                changeWorkspace(this);
                %getFromWorkspace(this);
            end
            
            % H81 Callback
            function btn_export_Callback(this, hObject, eventdata)
                handles = this.handles;
                % hObject    handle to btn_export (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                
                try
                    matRad_exportWidget;
                catch ME
                    handles = showError(this,'Could not export data.  Reason: ', ME);
                end
                this.handles = handles;
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
                    handles = showError(this,'Dose Import: Could not import data.  Reason: ', ME);
                    this.handles = handles;
                    return;
                end
                this.handles = handles;
                changeWorkspace(this);
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
                    for i = 1:length(RefVarNames)
                        if sum(ismember(AllVarNames,RefVarNames{i}))>0
                            evalin('base',['clear ', RefVarNames{i}]);
                        end
                    end
                    
                    %call the gui
                    h=matRad_importWidget;
                    uiwait(h.widgetHandle);
                    
                    %Check if we have the variables in the workspace
                    if evalin('base','exist(''cst'',''var'')') == 1 && evalin('base','exist(''ct'',''var'')') == 1
                        cst = evalin('base','cst');
                        ct = evalin('base','ct');
                        cst = generateCstTable(this,cst);
%                         handles.TableChanged = false;
%                         set(handles.popupTypeOfPlot,'Value',1);
                        
                        % compute HU values
                        if ~isfield(ct, 'cubeHU')
                            ct = matRad_electronDensitiesToHU(ct);
                            assignin('base','ct',ct);
                        end
%                         if ~isfield(ct, 'cubeHU')
%                             handles.cubeHUavailable = false;
%                         else
%                             handles.cubeHUavailable = true;
%                         end
                        
%                         % precompute contours
%                         cst = precomputeContours(this,ct,cst);
                        
                        assignin('base','ct',ct);
                        assignin('base','cst',cst);
                        
                        if evalin('base','exist(''pln'',''var'')')
                            assignin('base','pln',pln);
%                             setPln(handles);
%                         else
%                             getPlnFromGUI(handles);
%                             setPln(handles);
                        end
                    end
                    
                catch ME
                    handles = showError(this,'Binary Patient Import: Could not import data.  Reason: ', ME);
                   	this.handles = handles;
                    getFromWorkspace(this);
                    return;
                end
                
                this.handles = handles;
                changeWorkspace(this);
                %getFromWorkspace(this);
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
    
    
