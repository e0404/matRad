function matRad_modCst(cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad gui to change the optimization parameters in cst struct
% 
% call
%   matRad_modCst(cst)
%
% input
%   cst:    matRad cst struct
%   
% output
%   this function automatically updates the input argument in the workspace
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numOfAddedConstraints = 0;
numOfRows = 0;
    % create figure
    structWindow = figure('Name','matRad VOI/dose/penalty','NumberTitle','off',...
        'units','normalized','outerposition',[0 0 1 1],'ToolBar','figure','CloseRequestFcn',@CloseCallbackX);
    
    data               = guidata(gcf);
    data.inputname     = inputname(1);
    data.cst           = cst;
    data.structWindow  = structWindow;
    data.horViewOffset = 0.06;
    data.structDelete  = [];
    data.deleteCnt     = 1;
    data.AddIndex      = [];
    guidata(gcf,data);

    dataSetText = uicontrol('Style', 'text','String', 'Currently set parameters',...
        'FontSize', 12, 'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.35 0.95 0.3 0.03],'FontSize',14);
    set(dataSetText,'TooltipString',sprintf('configure your treatment plan')) ;
    
    dataSetText = uicontrol('Style', 'text', 'String', 'VOI','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.01+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('determines the volume of interest')) ;
     
    dataSetText = uicontrol('Style', 'text', 'String', 'VOI type','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.1+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('determines the tpye of the \n selected volume of interest')) ;
                 
    dataSetText = uicontrol('Style', 'text', 'String', 'Obj. func.',...
        'FontSize', 11, 'units', 'normalized', 'BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.3+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('determines the objective function for optimization')) ;
    
    dataSetText = uicontrol('Style', 'text', 'String', 'Overlap priority','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.2+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('determines the overlap priority of the volume of interest\n 0 = highest priority & inf = lowest priority')) ;
    
    dataSetText = uicontrol('Style', 'text', 'String', 'Penalty', 'FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.4+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('determines the corresponding penalty')) ;
                                       
    dataSetText = uicontrol('Style', 'text', 'String', 'Dose','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.5+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('prescripbed dose')) ;
                    
    dataSetText = uicontrol('Style', 'text', 'String', 'Exponent for EUD',...
        'FontSize', 11, 'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.60+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('Exponent in case of EUD objective function')) ;
    
    dataSetText = uicontrol('Style', 'text', 'String',' delete obj. func.',...
        'FontSize', 11, 'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.70+data.horViewOffset 0.9 0.075 0.03]);
    set(dataSetText,'TooltipString',sprintf('if enabled objective function will be used for optimization \n if disabled objective function will not be considered')) ;
    
    objectiveCounter = 1;
    
    % loop over all VOIs
    for i = 1:size(data.cst,1)
        % display only not ignored VOI
        if ~isequal(data.cst{i,3},'IGNORED')
            
            % loop over the number of objectives for the current VOI 
            k=1;
            while 1

                data.VOIText(objectiveCounter) = uicontrol('Style', 'text', 'String', data.cst{i,2},'FontSize', 10,...
                    'units', 'normalized','Position', [0.01+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02]);

                if isequal(data.cst{i,3}, 'TARGET')
                    BodyType = 1;
                else
                    BodyType = 2;
                end

                data.popupVOIType(objectiveCounter) = uicontrol('Style', 'popup', 'String', {'TARGET', 'OAR'},...
                    'FontSize', 10,'units', 'normalized','Position', [0.1+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                    'Callback', @popUpVOITypeCallback, 'Value', BodyType,'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter));
                
                if isempty(data.cst{i,6})
                    ObjFunc = 1;
                    data.cst{i,6}(k).parameter(1)=nan;
                    data.cst{i,6}(k).parameter(2)=nan;
                    data.cst{i,6}(k).exponent=nan;
                    data.cst{i,6}(k).type ='';
                elseif isequal(data.cst{i,6}(k).type, 'square underdosing')
                    ObjFunc = 2;
                elseif isequal(data.cst{i,6}(k).type, 'square overdosing')
                    ObjFunc = 3;
                elseif isequal(data.cst{i,6}(k).type, 'square deviation')
                    ObjFunc = 4;
                elseif isequal(data.cst{i,6}(k).type, 'mean')
                    ObjFunc = 5;
                elseif isequal(data.cst{i,6}(k).type, 'EUD')
                    ObjFunc = 6;
                end
                
                data.popupObjFunc(1,objectiveCounter) = uicontrol('Style', 'popup',...
                    'String',{'Please select ...','square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
                    'FontSize', 10, 'units', 'normalized','Position', [0.3+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],'Value', ObjFunc,...
                    'Callback', @popupObjFuncCallback, 'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter));
                
                data.editPriority(1,objectiveCounter) = uicontrol('Style', 'edit', 'String', data.cst{i,5}.Priority,...
                     'FontSize', 10, 'units', 'normalized','Position', [0.2+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                     'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter), 'Callback', @editPriorityCallback);    
                
                data.editPenalty(1,objectiveCounter) = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).parameter(1),...
                     'FontSize', 10, 'units', 'normalized','Position', [0.4+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                     'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter), 'Callback', @editPenaltyCallback);    
                
                if (~isequal(data.cst{i,6}(k).type, 'mean') && ~isequal(data.cst{i,6}(k).type, 'EUD')) 
                    
                    data.editDose(1,objectiveCounter) = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).parameter(2),...
                        'FontSize', 10, 'units', 'normalized','Position', [0.5+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                        'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter),'Callback', @editDoseCallback);
                    
                end
                 
                if isequal(data.cst{i,6}(k).type, 'EUD')
                    
                    data.editExponent(1,objectiveCounter) = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).exponent,...
                        'FontSize', 10, 'units', 'normalized','Position', [0.6+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                        'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter),'Callback', @editExponentCallback);
                    
                end
              
                data.btnDelete(1,objectiveCounter) = uicontrol('Style', 'pushbutton','String','delete',...
                    'FontSize', 10, 'units', 'normalized','Position',[0.7+data.horViewOffset 0.9-objectiveCounter*0.03 0.075 0.02],...
                    'Callback', @pushbtnDeleteCallback,'Tag', sprintf('%d,%d,%d', i,k,objectiveCounter)); 
                
                objectiveCounter = objectiveCounter+1;
                
                % necesarry to ensure behaviour of do-while loop
                  if k>=size(data.cst{i,6},2)
                    break
                  else
                      k=k+1;
                  end

            end
        end
        

    end
    numOfRows = numOfAddedConstraints;
    guidata(gcf, data);
    
    % accept and exit
    pushbuttonAccept = uicontrol('Style', 'pushbutton', 'String', 'Save & exit',...
        'FontSize',12, 'units', 'normalized', 'Position', [0.525 0.1 0.15 0.03],...
        'Callback', @pushbuttonAcceptCallback);
    set(pushbuttonAccept,'TooltipString',sprintf('Clicking this button will save the data and close the GUI')) ;
    
    % discard and exit
    pushbuttonExit = uicontrol('Style', 'pushbutton', 'String', 'Exit',...
        'FontSize',12, 'units', 'normalized', 'Position', [0.775 0.1 0.15 0.03],...
        'Callback', @CloseCallbackX);
    set(pushbuttonExit,'TooltipString',sprintf('Clicking this button will discard all changes and close the GUI')) ;
    
    % add a VOI or further objectives for a VOI
    pushbuttonAdd = uicontrol('Style', 'pushbutton', 'String', 'Add',...
        'FontSize', 12, 'units', 'normalized', 'Position', [0.375 0.1 0.1 0.03],...
        'Tag', sprintf('%d',objectiveCounter),'Callback', @pushbuttonAddCallback);
    set(pushbuttonAdd,'TooltipString',sprintf('Click here to add another objective function of your choice')) ;
 
uiwait

%% Definition Callback

    function popupObjFuncCallback(hObj, ~)
        tag = str2num(get(hObj,'Tag'));
        val = get(hObj, 'Value');
        before = data.cst{tag(1),6}(tag(2)).type;
        
        if (~isequal(before, 'mean') && ~isequal(before, 'EUD') && val == 5)
            
            set(data.editDose(tag(3)),'String','');
            set(data.editDose(tag(3)),'Visible','off');          
        
        elseif (~isequal(before, 'mean') && ~isequal(before,'EUD') && val ==6)
            
            set(data.editDose(tag(3)),'String','');
            set(data.editDose(tag(3)),'Visible','off'); 
            
            data.editExponent(tag(3)) = uicontrol('Style', 'edit',...
                'FontSize', 10, 'units', 'normalized','Position', [0.6+data.horViewOffset 0.9-tag(3)*0.03 0.075 0.02],...
                'Tag', sprintf('%d,%d,%d',tag(1),tag(2),tag(3)),'Callback', @editExponentCallback);
            
        elseif (isequal(before, 'EUD') && val ~= 5 && val ~= 6)
            
            data.editDose(tag(3)) = uicontrol('Style', 'edit', 'units', 'normalized',...
                'FontSize', 10, 'Position', [0.5+data.horViewOffset 0.9-tag(3)*0.03 0.075 0.02],...
                'Callback', @editDoseCallback, 'Tag', sprintf('%d,%d,%d',tag(1),tag(2),tag(3)));
            
            set(data.editExponent(tag(3)),'Visible','off'); 
            
        elseif (isequal(before, 'mean') && val ~= 5 && val ~= 6)
            
            set(data.editDose(tag(3)),'Visible','on');            
            
        elseif (isequal(before, 'EUD') && val == 5)
                        
            set(data.editExponent(tag(3)),'Visible','off'); 
            
        elseif (isequal(before, 'mean') && val == 6)
            
            data.editExponent(tag(3)) = uicontrol('Style', 'edit',...
                'FontSize', 10, 'units', 'normalized','Position', [0.6+data.horViewOffset 0.9-tag(3)*0.03 0.075 0.02],...
                'Tag', sprintf('%d,%d,%d',tag(1),tag(2),tag(3)),'Callback', @editExponentCallback);
                
        end
   
        guidata(gcf, data);

        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);            

        if val == 2
            data.cst{tag(1),6}(tag(2)).type = 'square underdosing';
        elseif val == 3
            data.cst{tag(1),6}(tag(2)).type = 'square overdosing';
        elseif val == 4
            data.cst{tag(1),6}(tag(2)).type = 'square deviation';
        elseif val == 5
            data.cst{tag(1),6}(tag(2)).type = 'mean';
        elseif val == 6
            data.cst{tag(1),6}(tag(2)).type = 'EUD';
        end
        guidata(gcf, data);
        
        function editExponentCallback(hObj, ~)
            tag = str2num(get(hObj,'Tag'));
            CheckValidity(str2num(get(hObj,'String')),data.editExponent(tag(3)));
        end
        
        function editDoseCallback(hObj, ~)
            tag = str2num(get(hObj, 'Tag'));
            CheckValidity(str2num(get(hObj,'String')),data.editDose(tag(3)));
        end
    
    end

    function editPriorityCallback(hObj, ~)
        
        data = guidata(gcf);
        
        tag = str2num(get(hObj,'Tag'));
        CheckValidity(str2num(get(hObj,'String')),data.editPriority(tag(3)));
                
        % change all priorities that belong to the same VOI
        list = get(data.VOIText(1,tag(3)),'String');
        value = get(data.VOIText(1,tag(3)),'Value');
        if iscell(list)
            CurrentVOI = list{value};
        else
            CurrentVOI = list;
        end
        
        CurrentPriority = str2num(get(hObj,'String'));

        for m = 1:size(data.VOIText,2)               
            list = get(data.VOIText(1,m),'String');
            value = get(data.VOIText(1,m),'Value');
            if iscell(list)
                LoopVOI = list{value};
            else
                LoopVOI = list;
            end

            if strcmp(CurrentVOI,LoopVOI)  
               set(data.editPriority(1,m),'String',num2str(CurrentPriority));
            end

        end
        
    end

    function editPenaltyCallback(hObj, ~)
        tag = str2num(get(hObj,'Tag'));
        CheckValidity(str2num(get(hObj,'String')),data.editPenalty(tag(3))); 
    end

    function editDoseCallback(hObj, ~)
        tag = str2num(get(hObj,'Tag'));
        CheckValidity(str2num(get(hObj,'String')),data.editDose(tag(3)));
    end

    function editExponentCallback(hObj, ~)
        tag = str2num(get(hObj,'Tag'));
        CheckValidity(str2num(get(hObj,'String')),data.editExponent(tag(3)));
    end

    function popUpVOITypeCallback(hObj, ~)          
        
        data = guidata(gcf);
        tag = str2num(get(hObj,'Tag')); 
        % change all priorities that belong to the same VOI
        list = get(data.VOIText(1,tag(3)),'String');
        value = get(data.VOIText(1,tag(3)),'Value');
        if iscell(list)
            CurrentVOI = list{value};
        else
            CurrentVOI = list;
        end
        % change all VOI Types that belong to the same VOI
        CurrVOIType = get(data.popupVOIType(1,tag(3)),'Value');
        for m = 1:size(data.VOIText,2) 
            if ~isnan(data.VOIText(1,m))
                 ListVOIText = get(data.VOIText(1,m),'String');
                 CurrentIndex = get(data.VOIText(1,m),'Value');
            else
                 ListVOIText=[];
                 CurrentIndex=[];
            end
            
            if iscell(ListVOIText)
                LoopVOI = ListVOIText{CurrentIndex};
            else
                LoopVOI = ListVOIText;
            end
            if strcmp(CurrentVOI,LoopVOI)                
                set(data.popupVOIType(1,m),'Value',CurrVOIType);
            end

        end
    end

    function pushbuttonAcceptCallback(~, ~) 
       
        % read data from gui, check for validity and write it into cst file
        FlagValidParameters = true;
        ObjFuncArray = get(data.popupObjFunc(1,find(sum(~isnan(data.popupObjFunc),1) > 0, 1 ,'first')),'String');
        NewCST = []; 
        Counter = 1;
  
        % loop over all rows in cst and gui and generate new cst
        for Cnt = 1:size(data.cst,1)
           
           CntObjF = 1;
           FlagFound = false;
           
           % loop over all rows in gui
           for j=1:length(data.editPenalty)   
                
               if ~isnan(data.VOIText(1,j))
                 ListVOIText = get(data.VOIText(1,j),'String');
                 CurrentIndex = get(data.VOIText(1,j),'Value');
               else
                   ListVOIText=[];
                   CurrentIndex=[];
               end
               
               if iscell(ListVOIText)
                  CurrTxtVOI = ListVOIText{CurrentIndex};
               else
                  CurrTxtVOI=ListVOIText;
               end
                
                    if strcmp(data.cst{Cnt,2},CurrTxtVOI) && strcmp(get(data.VOIText(1,j),'Enable'),'on')
                           
                           FlagFound = true;
                           NewCST{Counter,3} = data.cst{Cnt,2};
                           %get type of voi
                           Content = get(data.popupVOIType(1,j),'String');
                           NewCST{Counter,4} = Content{get(data.popupVOIType(1,j),'Value')};
                           %get objective function
                           NewCST{Counter,1}(CntObjF).type = ObjFuncArray{get(data.popupObjFunc(1,j),'Value'),1};
                           if strcmp(NewCST{Counter,1}(CntObjF).type,'Please select ...')
                               set(data.popupObjFunc(1,j),'BackgroundColor','r');
                               FlagValidParameters=false;
                           end

                           % get penalty
                           if isempty((get(data.editPenalty(1,j),'String'))) || isnan(str2num(get(data.editPenalty(1,j),'String')))
                                NewCST{Counter,1}(CntObjF).parameter=[];
                                set(data.editPenalty(1,j),'BackgroundColor','r');
                                FlagValidParameters=false;  
                           else
                                NewCST{Counter,1}(CntObjF).parameter(1,1)=str2num(get(data.editPenalty(1,j),'String'));
                           end
                           
                           
                           % get priority
                           NewCST{Counter,2}.Priority=str2num(get(data.editPriority(1,j),'String'));
                           if isempty(NewCST{Counter,2}.Priority) || isnan(NewCST{Counter,2}.Priority)
                                set(data.editPriority(1,j),'BackgroundColor','r');
                                FlagValidParameters=false;
                           end
                           
                           SelObjFunc = get(data.popupObjFunc(j),'Value');
                               
                           switch SelObjFunc

                               case 6
                                    NewCST{Counter,1}(CntObjF).exponent =str2num(get(data.editExponent(1,j),'String'));
                                    if isempty(NewCST{Counter,1}(CntObjF).exponent) || isnan(NewCST{Counter,1}(CntObjF).exponent)
                                        set(data.editExponent(1,j),'BackgroundColor','r');
                                        FlagValidParameters=false;
                                    end

                               case 5
                                    % do nothing in
                                    % this case
                                    
                               case {1,2,3,4} 

                                   if isempty((get(data.editDose(1,j),'String'))) || isnan(str2num(get(data.editDose(1,j),'String')))
                                        NewCST{Counter,1}(CntObjF).parameter=[];
                                        set(data.editDose(1,j),'BackgroundColor','r');
                                        FlagValidParameters=false;  
                                   else
                                        NewCST{Counter,1}(CntObjF).parameter(1,2)=str2num(get(data.editDose(1,j),'String'));
                                   end
                           end

                           CntObjF = CntObjF+1; 
                    end
                       
           end           
          
            if FlagFound == true
               Counter = Counter +1;
            end
           
        end

       % delete objectives in existing cst
       vIndex2Del = ones(size(data.cst,1),1);
       if FlagValidParameters
            for IdxCst = 1:size(data.cst,1)
                for ObjCnt = 1:size(data.cst{IdxCst,6},2)
                   for DelCnt = 1:size(data.structDelete,2)                  
                        if strcmp(data.cst(IdxCst,2),data.structDelete(DelCnt).VOI) ...
                                && strcmp(data.cst{IdxCst,6}(ObjCnt).type,data.structDelete(DelCnt).ObjFunc)                          
                                vIndex2Del(IdxCst)=0;
                        end
                   end
                end
            end
            
            Counter = 1 ;
            for IdxCst = 1:size(data.cst,1)
               if vIndex2Del(IdxCst)
                  tmpCst(Counter,:)=data.cst(IdxCst,:);
                  tmpCst{Counter,1}=Counter-1;
                  Counter = Counter + 1;
               end
            end
       end
       
       % write new cst in existing cst 
       if FlagValidParameters
            data.cst = [];
            data.cst = tmpCst;
            Counter = 1;
            for IdxCst = 1:size(NewCST,1)
                for IdxNewCst = 1:size(data.cst,1)
                    if  strcmp(data.cst(IdxNewCst,2),NewCST(IdxCst,3))
                        data.cst(Counter,6) = NewCST(IdxCst,1);
                        data.cst(Counter,3) = NewCST(IdxCst,4);
                        if isfield(NewCST{IdxCst,2},'Priority')
                            data.cst{Counter,5}.Priority = NewCST{IdxCst,2}.Priority;
                        else
                            data.cst{Counter,5}.Priority=nan;
                        end
                        Counter = Counter+1;
                    end
                    
                end
            end
            guidata(gcf, data);
            assignin('base',data.inputname,data.cst);
            CloseCallback();
        else
            warndlg('Empty fields detected');
        end
        
    end

    function pushbuttonAddCallback(hObj, ~)
        rowIdx = str2num(get(hObj,'Tag'))+numOfRows;
        
        tmp = data.AddIndex;
        tmp(end+1) = rowIdx;
        data.AddIndex = tmp;
        
        StringArray = cell(length(cst(:,2))+1,1);
        StringArray{1,1}='Select VOI..';
        StringArray(2:end,1)=cst(:,2);
        data.VOIText(1,rowIdx) = uicontrol('Style', 'popup', 'String', StringArray,'FontSize', 10,...
            'units', 'normalized', 'Position', [0.01+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
            'Callback', @popupVOICallback, 'Tag', sprintf('%d',rowIdx));
        
        function popupVOICallback(hObj, ~)
           
            VOI = get(hObj, 'Value')-1;
            rowIdx = str2num(get(hObj,'Tag'));
            
            % try to get VOI type of existing VOI
            list = get(data.VOIText(1,rowIdx),'String');
            value = get(data.VOIText(1,rowIdx),'Value');
            if iscell(list)
                CurrentVOI = list{value};
            else
                CurrentVOI = list;
            end
            ValueVOIType = 2;
            DefaultPriority = 1;
            % change all VOI Types that belong to the same VOI
            for m = 1:size(data.VOIText,2)-1               
                
                if ~isnan(data.VOIText(1,m))
                    list = get(data.VOIText(1,m),'String');
                    value = get(data.VOIText(1,m),'Value');
                else
                   list=[];
                   value=[];
                end
                if iscell(list)
                    LoopVOI = list{value};
                else
                    LoopVOI = list;
                end
                if strcmp(CurrentVOI,LoopVOI) && m<=size(data.popupVOIType,2)
                    ValueVOIType = get(data.popupVOIType(1,m),'Value');
                    DefaultPriority=str2num(get(data.editPriority(1,m),'String'));
                end

            end
            
            data.popupVOIType(1,rowIdx) = uicontrol('Style', 'popup', 'String', {'TARGET', 'OAR'},...
                    'Value', ValueVOIType,'FontSize', 10,'units', 'normalized', 'Position',[0.1+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                    'Callback', @popUpVOITypeCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
                
          
            data.popupObjFunc(1,rowIdx) = uicontrol('Style', 'popup',...
                'String', {'Please select ...','square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
                'FontSize',10,'units', 'normalized', 'Position', [0.3+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                'Callback', @popupObjFuncAddCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
                      
            data.editPriority(1,rowIdx) = uicontrol('Style', 'edit', 'FontSize', 10, ...
                'units', 'normalized', 'Position', [0.2+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                'String',DefaultPriority,...
                'Callback', @editPriorityCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
            
            data.editPenalty(1,rowIdx) = uicontrol('Style', 'edit', 'FontSize', 10, ...
                'units', 'normalized', 'Position', [0.4+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                'Callback', @editPenaltyCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
             
           data.editDose(1,rowIdx) = uicontrol('Style', 'edit', 'FontSize', 10,...
                 'units', 'normalized', 'Position', [0.5+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                 'Callback', @editDoseCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
                    
           data.editExponent(1,rowIdx) = uicontrol('Style', 'edit', 'FontSize', 10,...
                  'units', 'normalized', 'Position', [0.6+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                  'Callback', @editExponentCallback, 'Tag', sprintf('%d,%d,%d', VOI,0,rowIdx));
              
           data.btnDelete(1,rowIdx) = uicontrol('Style', 'pushbutton','String','delete',...
                    'FontSize', 10, 'units', 'normalized','Position',[0.7+data.horViewOffset 0.9-rowIdx*0.03 0.075 0.02],...
                    'Callback', @pushbtnDeleteCallback,'Tag', sprintf('%d,%d,%d', i,k,rowIdx)); 
           
            function popupObjFuncAddCallback(hObj, ~)
                tag = str2num(get(hObj,'Tag'));
                val = get(hObj, 'Value');
              
                if val == 2 || val == 3 || val == 4 

                    set(data.editDose(1,tag(3)),'Visible','on');
                    set(data.editExponent(1,tag(3)),'Visible','off');

                elseif val == 5
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'mean';
                    set(data.editDose(tag(3)),'String','');
                    set(data.editExponent(tag(3)),'String','');
                    set(data.editDose(1,tag(3)),'Visible','off');
                    set(data.editExponent(1,tag(3)),'Visible','off');                    
                    
                elseif val == 6
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'EUD';
                    set(data.editDose(tag(3)),'String','');
                    set(data.editDose(1,tag(3)),'Visible','off');
                    set(data.editExponent(1,tag(3)),'Visible','on');                         
                end
                
                if val > 1
                     set(hObj,'BackgroundColor',[0.94 0.94 0.94]);
                else
                     set(hObj,'BackgroundColor','r');
                end 
                
                guidata(gcf, data);             
            end
            
            guidata(gcf,data);
        end
       
        numOfAddedConstraints = numOfAddedConstraints + 1;
        numOfRows=numOfRows+1;
        guidata(gcf,data);
    end

    function pushbtnDeleteCallback(hObj, ~)
        
        tag = str2num(get(hObj,'Tag'));%Nummer des VOI, Nummer der Bedingung, Nummer der Zeile in Figure

        delete(data.btnDelete(tag(3)));
        data.btnDelete(tag(3)) = nan;
        
        if isfield(data,'editDose')
           if tag(3)<= length(data.editDose)
               delete(data.editDose(tag(3)))
               data.editDose(tag(3)) = nan;
           end
        end

        if isfield(data,'editExponent')
            if tag(3)<= length(data.editExponent)
                 delete(data.editExponent(tag(3)))
                 data.editExponent(tag(3))=nan;
            end
        end
               
        delete(data.editPenalty(tag(3)));
        data.editPenalty(tag(3))=nan;
        delete(data.editPriority(tag(3)));
        data.editPriority(tag(3))=nan;
        
        % remember VOI in order to delete from cst
        list = get(data.VOIText(1,tag(3)),'String');
        value = get(data.VOIText(1,tag(3)),'Value');
        if  iscell(list)
            CurrentVOI = list{value};
        else
            CurrentVOI = list;
        end
        
        
        % remember corresponding obj function in order to delete from cst
        list = get(data.popupObjFunc(1,tag(3)),'String');
        value = get(data.popupObjFunc(1,tag(3)),'Value');
        
        % only consider existing objective functions because added
        % contraints arent stored yet in cst
        if value>1 && isempty(data.AddIndex) || value>1 && sum(ismember(data.AddIndex,tag(3)))<1

            data.structDelete(data.deleteCnt).VOI=CurrentVOI;
            data.structDelete(data.deleteCnt).ObjFunc=list{value};
            data.deleteCnt = data.deleteCnt+1;
            
        end
   
        delete(data.VOIText(1,tag(3)));
        data.VOIText(tag(3))=nan;
        
        delete(data.popupObjFunc(tag(3)));
        data.popupObjFunc(tag(3))=nan;
        
        delete(data.popupVOIType(tag(3)));
        data.popupVOIType(tag(3))=nan;
        numOfAddedConstraints = numOfAddedConstraints-1;
        guidata(gcf,data);
    end

    function FlagValidity = CheckValidity(Val,hObj) 
      
       FlagValidity = true;
       
      if  isempty(Val)
               warndlg('Input not a number !');
               FlagValidity = false;        
      end
      
      if Val < 0 
          warndlg('Input has to be a positive number !');
          FlagValidity = false;  
      end
      
      if FlagValidity ==false
              set(hObj,'BackgroundColor','r');
        else
              set(hObj,'BackgroundColor',[0.94 0.94 0.94]);
      end
      
    end

    function CloseCallback(~, ~)
       
        selection = questdlg({'Are you done setting up objectives for optimization?','Changes have been saved.'},...
          'Closing Configuration Form',...
          'Yes','No','Yes'); 
        switch selection, 
          case 'Yes',
               delete(data.structWindow); 
               
          case 'No'
              return

        end
    end

    function CloseCallbackX(~, ~)
       
        selection = questdlg({'Really want to quit?','Changes will not be saved.'},...
          'Closing Configuration Form',...
          'Yes','No','Yes'); 
        switch selection, 
          case 'Yes',
               delete(data.structWindow); 
               
          case 'No'
              return

        end
    end

uiresume
close
 
end