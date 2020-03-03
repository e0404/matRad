classdef matRad_WorkflowWidget < matRad_Widget
    
    properties
    end
    
    methods
        function this = matRad_WorkflowWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[138.4 -7.38461538461539 273.4 59.5384615384615],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],... 'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','matRadGUI',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);
        end
       
        function changeWorkspace(obj)
            notify(obj, 'workspaceChanged');
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
        end
        
        methods
            
            % H74 Callback
            function btnLoadMat_Callback(this, hObject, event)
                % hObject    handle to btnLoadMat (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                
                handles = this.handles;
                
                [FileName, FilePath] = uigetfile('*.mat');
                if FileName == 0 % user pressed cancel --> do nothing.
                    return;
                end
                
                handles = resetGUI(this, hObject, event); ...resetGUI(hObject, handles, varargin)
                
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
                    %set(handles.legendTable,'String',{'no data loaded'});
                    %set(handles.popupDisplayOption,'String','no option available');
                    
                catch ME
                    handles = showError(this,'LoadMatFileFnc: Could not load *.mat file',ME);
                    
                    %    guidata(hObject,handles);
                    
                    %UpdateState(handles);
                    %UpdatePlot(handles);
                    
                    this. handles = handles;
                    return
                end
                
                try
                    generateCstTable(handles,cst);
                    handles.TableChanged = false;
                    set(handles.popupTypeOfPlot,'Value',1);
                    cst = matRad_computeVoiContoursWrapper(cst,ct);
                    
                    assignin('base','ct',ct);
                    assignin('base','cst',cst);
                    handles.State = 1;
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
                % if exist('stf','var') && exist('dij','var')
                %     handles.State = 2;
                % end
                
                if exist('resultGUI','var')
                    assignin('base','resultGUI',resultGUI);
                    % handles.State = 3;
                    % handles.SelectedDisplayOption ='physicalDose';
                end
                
                % recheck current workspace variables
                AllVarNames = evalin('base','who');
                handles.AllVarNames = AllVarNames;
                
                if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
                    handles = reloadGUI(hObject, handles, ct, cst);
                else
                    handles = reloadGUI(hObject, handles);
                end
                
                notify(this,'workspaceChanged');
                
                % guidata(hObject,handles);
                this.handles = handles;
                
                
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
                    
                    %pause(0.1);
                    %uiTable_CellEditCallback(hObject,[],handles);
                    %pause(0.3);
                    
                    %% get cst from table
                    %if ~getCstTable(handles)
                    %    return
                    %end
                    % read plan from gui and save it to workspace
                    % gets also IsoCenter from GUI if checkbox is not checked
                    getPlnFromGUI(handles);
                    
                    % get default iso center as center of gravity of all targets if not
                    % already defined
                    pln = evalin('base','pln');
                    
                    if length(pln.propStf.gantryAngles) ~= length(pln.propStf.couchAngles)
                        handles = showWarning(handles,'number of gantryAngles != number of couchAngles');
                    end
                    %%
                    if ~checkRadiationComposition(handles);
                        fileName = [pln.radiationMode '_' pln.machine];
                        handles = showError(handles,errordlg(['Could not find the following machine file: ' fileName ]));
                        %guidata(hObject,handles);
                        this.handles = handles;
                        return;
                    end
                    
                    %% check if isocenter is already set
                    if ~isfield(pln.propStf,'isoCenter')
                        handles = showWarning(handles,'no iso center set - using center of gravity based on structures defined as TARGET');
                        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
                        assignin('base','pln',pln);
                    elseif ~get(handles.checkIsoCenter,'Value')
                        if ~strcmp(get(handles.editIsoCenter,'String'),'multiple isoCenter')
                            pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*str2num(get(handles.editIsoCenter,'String'));
                        end
                    end
                    
                catch ME
                    handles = showError(handles,'CalcDoseCallback: Error in preprocessing!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    % guidata(hObject,handles);
                    this.handles = handles;
                    return;
                end
                
                % generate steering file
                try
                    currPln = evalin('base','pln');
                    % if we run 3d conf opt -> hijack runDao to trigger computation of
                    % connected bixels
                    if strcmp(pln.radiationMode,'photons') && get(handles.radiobutton3Dconf,'Value')
                        currpln.propOpt.runDAO = true;
                    end
                    stf = matRad_generateStf(evalin('base','ct'),...
                        evalin('base','cst'),...
                        currPln);
                    assignin('base','stf',stf);
                catch ME
                    handles = showError(handles,'CalcDoseCallback: Error in steering file generation!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    % guidata(hObject,handles);
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
                    handles.State = 2;
                    handles.TableChanged = false;
                    UpdateState(handles);
                    UpdatePlot(handles);
                    % guidata(hObject,handles);
                    this.handles = handles;
                    
                catch ME
                    handles = showError(handles,'CalcDoseCallback: Error in dose calculation!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    % guidata(hObject,handles);
                    this.handles = handles;
                    return;
                end
                
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                
                % guidata(hObject,handles);
                this.handles = handles;
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
                    btnTableSave_Callback([],[],handles); %We don't need it?
                    
                    
                    % if a critical change to the cst has been made which affects the dij matrix
                    if handles.DijCalcWarning == true
                        
                        choice = questdlg('Overlap priorites of OAR constraints have been edited, a new OAR VOI was added or a critical row constraint was deleted. A new Dij calculation might be necessary.', ...
                            'Title','Cancel','Calculate Dij then Optimize','Optimze directly','Optimze directly');
                        
                        switch choice
                            case 'Cancel'
                                set(Figures, 'pointer', 'arrow');
                                set(InterfaceObj,'Enable','on');
                                % guidata(hObject,handles);
                                this.handles = handles;
                                
                                return;
                            case 'Calculate dij again and optimize'
                                handles.DijCalcWarning = false;
                                btnCalcDose_Callback(hObject, eventdata, handles)
                            case 'Optimze directly'
                                handles.DijCalcWarning = false;
                        end
                    end
                    
                    pln = evalin('base','pln');
                    ct  = evalin('base','ct');
                    
                    % optimize
                    if get(handles.radiobutton3Dconf,'Value') && strcmp(handles.Modalities{get(handles.popupRadMode,'Value')},'photons')
                        % conformal plan if photons and 3d conformal
                        if ~matRad_checkForConnectedBixelRows(evalin('base','stf'))
                            error('disconnetced dose influence data in BEV - run dose calculation again with consistent settings');
                        end
                        [resultGUIcurrentRun,usedOptimizer] = matRad_fluenceOptimization(matRad_collapseDij(evalin('base','dij')),evalin('base','cst'),pln);
                        resultGUIcurrentRun.w = resultGUIcurrentRun.w * ones(evalin('base','dij.totalNumOfBixels'),1);
                        resultGUIcurrentRun.wUnsequenced = resultGUIcurrentRun.w;
                    else
                        if pln.propOpt.runDAO
                            if ~matRad_checkForConnectedBixelRows(evalin('base','stf'))
                                error('disconnetced dose influence data in BEV - run dose calculation again with consistent settings');
                            end
                        end
                        
                        [resultGUIcurrentRun,usedOptimizer] = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),pln);
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
                    
                    % set some values
                    if handles.plane == 1
                        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,2)/ct.resolution.x));
                    elseif handles.plane == 2
                        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,1)/ct.resolution.y));
                    elseif handles.plane == 3
                        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,3)/ct.resolution.z));
                    end
                    
                    handles.State = 3;
                    handles.SelectedDisplayOptionIdx = 1;
                    if strcmp(pln.radiationMode,'carbon') || (strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD'))
                        handles.SelectedDisplayOption = 'RBExDose';
                    else
                        handles.SelectedDisplayOption = 'physicalDose';
                    end
                    handles.selectedBeam = 1;
                    % check IPOPT status and return message for GUI user if no DAO or
                    % particles
                    if ~pln.propOpt.runDAO || ~strcmp(pln.radiationMode,'photons')
                        CheckOptimizerStatus(usedOptimizer,'Fluence')
                    end
                    
                catch ME
                    handles = showError(handles,'OptimizeCallback: Could not optimize!',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    % guidata(hObject,handles);
                    this.handles = handles;
                    return;
                end
                
                % perform sequencing and DAO
                try
                    
                    %% sequencing
                    if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
                        %   resultGUI = matRad_xiaLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
                        %       ,get(handles.editSequencingLevel,'Value'));
                        %   resultGUI = matRad_engelLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
                        %       ,str2double(get(handles.editSequencingLevel,'String')));
                        resultGUI = matRad_siochiLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
                            ,str2double(get(handles.editSequencingLevel,'String')));
                        
                        assignin('base','resultGUI',resultGUI);
                    end
                    
                catch ME
                    handles = showError(handles,'OptimizeCallback: Could not perform sequencing',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    %guidata(hObject,handles);
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
                        CheckOptimizerStatus(usedOptimizer,'DAO');
                    end
                    
                    if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
                        matRad_visApertureInfo(resultGUI.apertureInfo);
                    end
                    
                catch ME
                    handles = showError(handles,'OptimizeCallback: Could not perform direct aperture optimization',ME);
                    % change state from busy to normal
                    set(Figures, 'pointer', 'arrow');
                    set(InterfaceObj,'Enable','on');
                    % guidata(hObject,handles);
                    this.handles = handles;
                    return;
                end
                
                % change state from busy to normal
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                handles.dispWindow{3,1}  = [];   % reset dose ranges
                handles.dispWindow{3,2}  = [];   % reset min max dose values
                handles.rememberCurrAxes = false;
                handles.IsoDose.Levels   = 0;  % ensure to use default iso dose line spacing
                handles.cBarChanged      = true;
                
                % guidata(hObject,handles);
                this.handles = handles;
                handles = updateIsoDoseLineCache(handles);
                UpdateState(handles);
                UpdatePlot(handles);
                handles.rememberCurrAxes = true;
                % guidata(hObject,handles);
                this.handles = handles;
            end
            
            % H77 Callback
            function btnLoadDicom_Callback(this, hObject, event)
                % hObject    handle to btnLoadDicom (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                handles = this.handles;
                try
                    % delete existing workspace - parse variables from base workspace
                    set(handles.popupDisplayOption,'String','no option available');
                    AllVarNames = evalin('base','who');
                    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
                    for i = 1:length(RefVarNames)
                        if sum(ismember(AllVarNames,RefVarNames{i}))>0
                            evalin('base',['clear ', RefVarNames{i}]);
                        end
                    end
                   % handles.State = 0;
                    matRad_importDicomWidget; % matRad_importDicomGUI;
                    
                catch
                    handles = showError(handles,'DicomImport: Could not import data');
                end
               % UpdateState(handles);
                %guidata(hObject,handles);
                this.handles = handles;
            end
            
            % H78 Callback - button: refresh
            function btnRefresh_Callback(this, hObject, event)
                handles = this.handles;
                handles = resetGUI(this, hObject, event); ...resetGUI(hObject, handles, varargin)
                %% parse variables from base workspace
                AllVarNames = evalin('base','who');
                handles.AllVarNames = AllVarNames;
                try
                    if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
                        ct  = evalin('base','ct');
                        cst = evalin('base','cst');
                        %cst = setCstTable(handles,cst);
                        generateCstTable(handles,cst);
                      %  handles.State = 1;
                        cst = matRad_computeVoiContoursWrapper(cst,ct);
                        assignin('base','cst',cst);
                    elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
                      %  handles = showError(handles,'GUI OpeningFunc: could not find cst file');
                    elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
                     %   handles = showError(handles,'GUI OpeningFunc: could not find ct file');
                    end
                catch
                  % handles = showError(handles,'GUI OpeningFunc: Could not load ct and cst file');
                end
                
                if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
                   handles = initialize(this);...  handles = reloadGUI(hObject, handles, ct, cst);
                else
                   handles = initialize(this); ... handles = reloadGUI(hObject, handles);
                end
                this.handles = handles;
            end
            
%             function btnRefresh_Callback(this, hObject, event)
%                 handles = this.handles;
%                % handles = resetGUI(hObject, handles);
%                 
%                 %% parse variables from base workspace
%                 AllVarNames = evalin('base','who');
%                 handles.AllVarNames = AllVarNames;
%                 try
%                     if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
%                         ct  = evalin('base','ct');
%                         cst = evalin('base','cst');
%                         %cst = setCstTable(handles,cst);
%                         generateCstTable(handles,cst);
%                       %  handles.State = 1;
%                         cst = matRad_computeVoiContoursWrapper(cst,ct);
%                         assignin('base','cst',cst);
%                     elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
%                         handles = showError(handles,'GUI OpeningFunc: could not find cst file');
%                     elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
%                      %   handles = showError(handles,'GUI OpeningFunc: could not find ct file');
%                     end
%                 catch
%                   % handles = showError(handles,'GUI OpeningFunc: Could not load ct and cst file');
%                 end
%                 
%                 if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
%                    handles = initialize(this);...  handles = reloadGUI(hObject, handles, ct, cst);
%                 else
%                    handles = initialize(this); ... handles = reloadGUI(hObject, handles);
%                 end
%                 %guidata(hObject, handles);
%                 this.handles = handles;
%             end
            
            % H79 Callback
            function pushbutton_recalc_Callback(this, hObject, eventdata)
                
                handles = this.handles;
                % recalculation only makes sense if ...
                if evalin('base','exist(''pln'',''var'')') && ...
                        evalin('base','exist(''stf'',''var'')') && ...
                        evalin('base','exist(''ct'',''var'')') && ...
                        evalin('base','exist(''cst'',''var'')') && ...
                        evalin('base','exist(''resultGUI'',''var'')')
                    
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
                        
                        % get weights of the selected cube
                        Content = get(handles.popupDisplayOption,'String');
                        SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
                        Suffix = strsplit(SelectedCube,'_');
                        if length(Suffix)>1
                            Suffix = ['_' Suffix{2}];
                        else
                            Suffix = '';
                        end
                        
                        if sum([stf.totalNumOfBixels]) ~= length(resultGUI.(['w' Suffix]))
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
                        resultGUIreCalc = matRad_calcCubes(resultGUI.(['w' Suffix]),dij,cst);
                        
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
                        
                        handles.State = 3;
                        
                        % show physicalDose of newly computed state
                        handles.SelectedDisplayOption = 'physicalDose';
                        set(handles.popupDisplayOption,'Value',find(strcmp('physicalDose',Content)));
                        
                        % change state from busy to normal
                        set(Figures, 'pointer', 'arrow');
                        set(InterfaceObj,'Enable','on');
                        
                        handles.cBarChanged = true;
                        
                        handles = updateIsoDoseLineCache(handles);
                        
                        UpdateState(handles);
                        
                        handles.rememberCurrAxes = false;
                        UpdatePlot(handles);
                        handles.rememberCurrAxes = true;
                        
                        % guidata(hObject,handles);
                        this.handles = handles;
                        
                    catch ME
                        handles = showError(handles,'CalcDoseCallback: Error in dose recalculation!',ME);
                        
                        % change state from busy to normal
                        set(Figures, 'pointer', 'arrow');
                        set(InterfaceObj,'Enable','on');
                        
                        % guidata(hObject,handles);
                        this.handles = handles;
                        return;
                        
                    end
                    
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
                        'Callback', @(hpb,eventdata)SaveResultToGUI(hpb,eventdata,guidata(hpb)));
                end
                
                uiwait(figDialog);
                % guidata(hObject, handles);
                this.handles = handles;
              %  UpdateState(handles)
                UpdatePlot(handles)
            end
            
            % H81 Callback
            function btn_export_Callback(this, hObject, eventdata)
                handles = this.handles;
                % hObject    handle to btn_export (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                
                try
                    % matRad_exportGUI;
                    matRad_exportWidget;
                catch
                    handles = showError(handles,'Could not export data');
                end
               % UpdateState(handles);
                % guidata(hObject,handles);
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
                btnRefresh_Callback(this, hObject, eventdata)
                this. handles = handles;
            end
            
             % H83 Callback
            function pushbutton_importFromBinary_Callback(this, hObject, eventdata)
                % hObject    handle to pushbutton_importFromBinary (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                handles = this.handles;
                
                try
                    % delete existing workspace - parse variables from base workspace
                    set(handles.popupDisplayOption,'String','no option available');
                    AllVarNames = evalin('base','who');
                    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
                    for i = 1:length(RefVarNames)
                        if sum(ismember(AllVarNames,RefVarNames{i}))>0
                            evalin('base',['clear ', RefVarNames{i}]);
                        end
                    end
                    handles.State = 0;
                    
                    %call the gui
                    uiwait(matRad_importGUI);
                    
                    %Check if we have the variables in the workspace
                    if evalin('base','exist(''cst'',''var'')') == 1 && evalin('base','exist(''ct'',''var'')') == 1
                        cst = evalin('base','cst');
                        ct = evalin('base','ct');
                        %setCstTable(handles,cst);
                        generateCstTable(hanles,cst);
                        handles.TableChanged = false;
                        set(handles.popupTypeOfPlot,'Value',1);
                        
                        % compute HU values
                        if ~isfield(ct, 'cubeHU')
                            ct = matRad_electronDensitiesToHU(ct);
                            assignin('base','ct',ct);
                        end
                        if ~isfield(ct, 'cubeHU')
                            handles.cubeHUavailable = false;
                        else
                            handles.cubeHUavailable = true;
                        end
                        
                        % precompute contours
                        cst = precomputeContours(ct,cst);
                        
                        assignin('base','ct',ct);
                        assignin('base','cst',cst);
                        
                        if evalin('base','exist(''pln'',''var'')')
                            assignin('base','pln',pln);
                            setPln(handles);
                        else
                            getPlnFromGUI(handles);
                            setPln(handles);
                        end
                       % handles.State = 1;
                    end
                    
                    % set slice slider
                    handles.plane = get(handles.popupPlane,'value');
                    %if handles.State >0
                        set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(handles.plane),...
                            'Value',round(ct.cubeDim(handles.plane)/2),...
                            'SliderStep',[1/(ct.cubeDim(handles.plane)-1) 1/(ct.cubeDim(handles.plane)-1)]);
                  %  end
                    
                   % if handles.State > 0
                        % define context menu for structures
                        for i = 1:size(cst,1)
                            if cst{i,5}.Visible
                                handles.VOIPlotFlag(i) = true;
                            else
                                handles.VOIPlotFlag(i) = false;
                            end
                        end
                    % end
                    
                    handles.dispWindow = cell(3,2);
                    handles.cBarChanged = true;
                    
                  %  UpdateState(handles);
                    handles.rememberCurrAxes = false;
                    UpdatePlot(handles);
                    handles.rememberCurrAxes = true;
                catch
                    handles = showError(handles,'Binary Patient Import: Could not import data');
                   % UpdateState(handles);
                end
                
                % guidata(hObject,handles);
                this.handles = handles;
                
            end
            
            function handles = resetGUI(this, hObject, eventdata) %(hObject, handles, varargin)
                
                handles = this.handles;
                
                [env, versionString] = matRad_getEnvironment();
                
                if strcmp(env,'MATLAB')
                    %OpenGL only works for pc (maybe also for linux?)
                    if ispc
                        opengl software
                    end
                end
                
            end
            
            
            function handles = reloadGUI(this,hObject, ct, cst)
                
                handles = this.handles;
                AllVarNames = handles.AllVarNames;
                
                if ismember('ct',AllVarNames)
                    % compute HU values
                    if ~isfield(ct, 'cubeHU')
                        ct = matRad_electronDensitiesToHU(ct);
                        assignin('base','ct',ct);
                    end
                    if ~isfield(ct, 'cubeHU')
                        handles.cubeHUavailable = false;
                    else
                        handles.cubeHUavailable = true;
                    end
                end
                
                this.handles = handles;
            end
        end
    end
    
    
