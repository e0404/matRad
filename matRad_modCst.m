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


global l

l = 0;

if nargin > 0
    
    % create figure
    structWindow = figure('Name','matRad VOI/dose/penalty','NumberTitle','off',...
        'units','normalized','outerposition',[0 0 1 1],'ToolBar','figure');
    
    data = guidata(gcf);
    data.inputname = inputname(1);
    data.cst = cst;
    guidata(gcf,data);
    
    dataSetText = uicontrol('Style', 'text','String', 'Currently set parameters',...
        'FontSize', 12, 'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.35 0.95 0.3 0.03]);
                  
    dataSetText = uicontrol('Style', 'text', 'String', 'VOI','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.03 0.9 0.075 0.03]);
     
    dataSetText = uicontrol('Style', 'text', 'String', 'VOI type','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.15 0.9 0.075 0.03]);
                 
    dataSetText = uicontrol('Style', 'text', 'String', 'Obj. func.',...
        'FontSize', 11, 'units', 'normalized', 'BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.27 0.9 0.075 0.03]);
                                       
    dataSetText = uicontrol('Style', 'text', 'String', 'Penalty', 'FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.39 0.9 0.075 0.03]);
                                       
    dataSetText = uicontrol('Style', 'text', 'String', 'Dose','FontSize', 11,...
        'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.51 0.9 0.075 0.03]);                   
                    
    dataSetText = uicontrol('Style', 'text', 'String', 'Exponent for EUD',...
        'FontSize', 11, 'units', 'normalized','BackgroundColor', [0.8 0.8 0.8],...
        'Position', [0.63 0.9 0.075 0.03]);
    
    objectiveCounter = 1;
    
    % loop over all VOIs
    for i = 1:size(data.cst,1)
        % display only not ignored VOI
        if ~isequal(data.cst{i,3},'IGNORED')
            
            dataSetText = uicontrol('Style', 'text', 'String', data.cst{i,2},'FontSize', 10,...
                'units', 'normalized','Position', [0.03 0.9-objectiveCounter*0.03 0.075 0.02]);
            
            if isequal(data.cst{i,3}, 'TARGET')
                BodyType = 1;
            else
                BodyType = 2;
            end
            
            popupBodyType = uicontrol('Style', 'popup', 'String', {'TARGET', 'OAR'},...
                'FontSize', 10,'units', 'normalized','Position', [0.15 0.9-objectiveCounter*0.03 0.075 0.02],...
                'Callback', @popupBodyTypeCallback, 'Value', BodyType);
            
            % loop over the number of objectives for the current VOI
            for k = 1:size(data.cst{i,6},2)
                
                if isequal(data.cst{i,6}(k).type, 'square underdosing')
                    ObjFunc = 1;
                elseif isequal(data.cst{i,6}(k).type, 'square overdosing')
                    ObjFunc = 2;
                elseif isequal(data.cst{i,6}(k).type, 'square deviation')
                    ObjFunc = 3;
                elseif isequal(data.cst{i,6}(k).type, 'mean')
                    ObjFunc = 4;
                elseif isequal(data.cst{i,6}(k).type, 'EUD')
                    ObjFunc = 5;
                end
                
                popupObjFunc(i) = uicontrol('Style', 'popup',...
                    'String',{'square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
                    'FontSize', 10, 'units', 'normalized','Position', [0.27 0.9-objectiveCounter*0.03 0.075 0.02],'Value', ObjFunc,...
                    'Callback', @popupObjFuncCallback, 'Tag', sprintf('%d,%d,%d',i,k,objectiveCounter));
                
                editPenalty = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).parameter(1),...
                     'FontSize', 10, 'units', 'normalized','Position', [0.39 0.9-objectiveCounter*0.03 0.075 0.02],...
                     'Tag', sprintf('%d,%d',i,k), 'Callback', @editPenaltyCallback);    
                
                if (~isequal(data.cst{i,6}(k).type, 'mean') && ~isequal(data.cst{i,6}(k).type, 'EUD')) 
                    
                    editDose = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).parameter(2),...
                        'FontSize', 10, 'units', 'normalized','Position', [0.51 0.9-objectiveCounter*0.03 0.075 0.02],...
                        'Tag', sprintf('%d,%d',i,k),'Callback', @editDoseCallback);
                    
                end
                 
                if isequal(data.cst{i,6}(k).type, 'EUD')
                    
                    editExponent = uicontrol('Style', 'edit', 'String', data.cst{i,6}(k).exponent,...
                        'FontSize', 10, 'units', 'normalized','Position', [0.63 0.9-objectiveCounter*0.03 0.075 0.02],...
                        'Tag', sprintf('%d,%d',i,k),'Callback', @editExponentCallback);
                    
                end
                
                pushbuttonDelete = uicontrol('Style', 'pushbutton', 'String', 'Delete',...
                    'FontSize', 10, 'units', 'normalized', 'Position', [0.75 0.9-objectiveCounter*0.03 0.075 0.02],...
                    'Callback', @pushbuttonDeleteCallback, 'Tag', sprintf('%d,%d,%d', i,k,objectiveCounter));
                
                objectiveCounter = objectiveCounter+1;
                
            end
        end
                      
    end
    
    % start matRad
    pushbuttonAccept = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
        'FontSize',12, 'units', 'normalized', 'Position', [0.525 0.1 0.1 0.03],...
        'Callback', @pushbuttonAcceptCallback);
    
    % add a VOI or further objectives for a VOI
    pushbuttonAdd = uicontrol('Style', 'pushbutton', 'String', 'Add',...
        'FontSize', 12, 'units', 'normalized', 'Position', [0.375 0.1 0.1 0.03],...
        'Tag', sprintf('%d',objectiveCounter),'Callback', @pushbuttonAddCallback);

end

%% Definition Callback

    function popupObjFuncCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'));
        val = get(hObj, 'Value');
        before = data.cst{tag(1),6}(tag(2)).type;
        
        if (~isequal(before, 'mean') && ~isequal(before, 'EUD') && val == 4)
            
            dose = uicontrol('Style', 'text', 'String', {}, 'units', 'normalized',...
                'Position', [0.51 0.9-tag(3)*0.03 0.075 0.02], 'BackgroundColor', [0.8 0.8 0.8]);
            
            % delete previous dose parameter
            data.cst{tag(1),6}(tag(2)).parameter(2)=[];
        
        elseif (~isequal(before, 'mean') && ~isequal(before,'EUD') && val ==5)
            
            dose = uicontrol('Style', 'text', 'String', {}, 'units', 'normalized',...
                'Position', [0.51 0.9-tag(3)*0.03 0.075 0.02], 'BackgroundColor', [0.8 0.8 0.8]);
            
            % delete previous dose parameter
            data.cst{tag(1),6}(tag(2)).parameter(2) = [];
            
            editExponent = uicontrol('Style', 'edit',...
                'FontSize', 10, 'units', 'normalized','Position', [0.63 0.9-tag(3)*0.03 0.075 0.02],...
                'Tag', sprintf('%d,%d',tag(1),tag(2)),'Callback', @editExponentCallback);
            
        elseif (isequal(before, 'EUD') && val ~= 4 && val ~= 5)
            
            editDose = uicontrol('Style', 'edit', 'units', 'normalized',...
                'FontSize', 10, 'Position', [0.51 0.9-tag(3)*0.03 0.075 0.02],...
                'Callback', @editDoseCallback, 'Tag', sprintf('%d,%d', tag(1), tag(2)));
            
            exponent = uicontrol('Style', 'text', 'String', {}, 'units', 'normalized',...
                'Position', [0.63 0.9-tag(3)*0.03 0.075 0.02], 'BackgroundColor', [0.8 0.8 0.8]);
            
        elseif (isequal(before, 'mean') && val ~= 4 && val ~= 5)
            
            editDose = uicontrol('Style', 'edit', 'units', 'normalized',...
                'FontSize', 10, 'Position', [0.51 0.9-tag(3)*0.03 0.075 0.02],...
                'Callback', @editDoseCallback, 'Tag', sprintf('%d,%d', tag(1), tag(2)));
            
            
        elseif (isequal(before, 'EUD') && val == 4)
                        
            exponent = uicontrol('Style', 'text', 'String', {}, 'units', 'normalized',...
                'Position', [0.63 0.9-tag(3)*0.03 0.075 0.02], 'BackgroundColor', [0.8 0.8 0.8]);
            
            %vorheriger Exponent wird gelöscht
            data.cst{tag(1),6}(tag(2)).exponent = [];
            
        elseif (isequal(before, 'mean') && val == 5)
            
            editExponent = uicontrol('Style', 'edit',...
                'FontSize', 10, 'units', 'normalized','Position', [0.63 0.9-tag(3)*0.03 0.075 0.02],...
                'Tag', sprintf('%d,%d',tag(1),tag(2)),'Callback', @editExponentCallback);
                
        end
        guidata(gcf, data);

        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);            

        if val == 1
            data.cst{tag(1),6}(tag(2)).type = 'square underdosing';
            data.cst{tag(1),6}(tag(2)).exponent = [];
        elseif val == 2
            data.cst{tag(1),6}(tag(2)).type = 'square overdosing';
            data.cst{tag(1),6}(tag(2)).exponent = [];
        elseif val == 3
            data.cst{tag(1),6}(tag(2)).type = 'square deviation';
            data.cst{tag(1),6}(tag(2)).exponent = [];
        elseif val == 4
            data.cst{tag(1),6}(tag(2)).type = 'mean';
            data.cst{tag(1),6}(tag(2)).exponent = [];
        elseif val == 5
            data.cst{tag(1),6}(tag(2)).type = 'EUD';
        end
        guidata(gcf, data);
        
        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);

        function editExponentCallback(hObj, event)
            tag = str2num(get(hObj,'Tag'));
            data.cst{tag(1),6}(tag(2)).exponent = str2num(get(hObj, 'String'));
            guidata(gcf, data);
            % overwrite cst struct directly in workspace
            assignin('base',data.inputname,data.cst);
        end
        
        function editDoseCallback(hObj, event)
            tag = str2num(get(hObj, 'Tag'));
            data.cst{tag(1),6}(tag(2)).parameter(2) = str2num(get(hObj,'String'));
            guidata(gcf, data);
            % overwrite cst struct directly in workspace
            assignin('base',data.inputname,data.cst);
        end
    
    end

    function editPenaltyCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'));
        data.cst{tag(1),6}(tag(2)).parameter(1) = str2num(get(hObj, 'String'));
        guidata(gcf,data);
        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);
    end

    function editDoseCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'));
        data.cst{tag(1),6}(tag(2)).parameter(2) = str2num(get(hObj, 'String'));
        guidata(gcf,data);
        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);
    end

    function editExponentCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'));
        data.cst{tag(1),6}(tag(2)).exponent = str2num(get(hObj, 'String'));
        guidata(gcf, data);
        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);
    end

    function pushbuttonAcceptCallback(hObj, event)
        close(structWindow);
        
    end

    function pushbuttonAddCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'))+l; %erste freie zeile
        popupVOI = uicontrol('Style', 'popup', 'String', cst(:,2),'FontSize', 10,...
            'units', 'normalized', 'Position', [0.03 0.9-tag*0.03 0.075 0.02],...
            'Callback', @popupVOICallback, 'Tag', sprintf('%d',tag));
        
        function popupVOICallback(hObj, event)
            
            VOI = get(hObj, 'Value');
            tag = str2num(get(hObj,'Tag'));
            
            if isequal(data.cst{VOI,3}, 'IGNORED')
                popupBodyType = uicontrol('Style', 'popup', 'String', {'TARGET', 'OAR'},...
                    'FontSize', 10,'units', 'normalized', 'Position',[0.15 0.9-tag*0.03 0.075 0.02],...
                    'Callback', @popupBodyTypeCallback, 'Tag', sprintf('%d', VOI));
                
            else
                dataSetText = uicontrol('Style', 'text', 'String', data.cst{VOI,3},...
                    'FontSize', 10, 'units', 'normalized', 'Position',[0.15 0.9-tag*0.03 0.075 0.02]);
            
            end
            
            popupObjFuncAdd = uicontrol('Style', 'popup',...
                'String', {'square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
                'FontSize',10,'units', 'normalized', 'Position', [0.27 0.9-tag*0.03 0.075 0.02],...
                'Callback', @popupObjFuncAddCallback, 'Tag', sprintf('%d, %d', VOI, tag));
            
            editPenaltyAdd = uicontrol('Style', 'edit', 'FontSize', 10, ...
                'units', 'normalized', 'Position', [0.39 0.9-tag*0.03 0.075 0.02],...
                'Callback', @editPenaltyAddCallback, 'Tag', sprintf('%d', VOI));
            
            function popupBodyTypeCallback(hObj, event)
            
                tag = str2num(get(hObj,'Tag')); %Zeilennummer des ausgewählten VOIs
                val = num2str(get(hObj, 'Value'));
                if val == '1',                
                    data.cst{tag,3} = 'TARGET';   
                else 
                    data.cst{tag,3} = 'OAR';
                end
                guidata(gcf, data);
                % overwrite cst struct directly in workspace
                assignin('base',data.inputname,data.cst);
            end
            
            function popupObjFuncAddCallback(hObj, event)
                tag = str2num(get(hObj,'Tag')); % number of voi
                val = get(hObj, 'Value');
                
                if val == 1
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'square underdosing';
                
                    editDoseAdd = uicontrol('Style', 'edit', 'FontSize', 10,...
                        'units', 'normalized', 'Position', [0.51 0.9-tag(2)*0.03 0.075 0.02],...
                        'Callback', @editDoseAddCallback, 'Tag', sprintf('%d', VOI));

                elseif val == 2
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'square overdosing';
                    
                    editDoseAdd = uicontrol('Style', 'edit', 'FontSize', 10,...
                        'units', 'normalized', 'Position', [0.51 0.9-tag(2)*0.03 0.075 0.02],...
                        'Callback', @editDoseAddCallback, 'Tag', sprintf('%d', VOI));
                    
                elseif val == 3
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'square deviation';
                    
                    editDoseAdd = uicontrol('Style', 'edit', 'FontSize', 10,...
                        'units', 'normalized', 'Position', [0.51 0.9-tag(2)*0.03 0.075 0.02],...
                        'Callback', @editDoseAddCallback, 'Tag', sprintf('%d', VOI));
                    
                elseif val == 4
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'mean';
                elseif val == 5
                    data.cst{tag(1),6}(size(data.cst{tag(1),6},2)+1).type = 'EUD';
                    
                    editExponentAdd = uicontrol('Style', 'edit', 'FontSize', 10,...
                        'units', 'normalized', 'Position', [0.63 0.9-tag(2)*0.03 0.075 0.02],...
                        'Callback', @editExponentAddCallback, 'Tag', sprintf('%d', VOI));
                end
                
                guidata(gcf, data);
    
                % overwrite cst struct directly in workspace
                assignin('base',data.inputname,data.cst);
                
                function editDoseAddCallback(hObj, event)
                    tag= str2num(get(hObj,'Tag')); %Nummer des VOIs
                    data.cst{tag,6}(size(data.cst{tag,6},2)).parameter(2) = str2num(get(hObj, 'String'));
                    guidata(gcf, data);
                    % overwrite cst struct directly in workspace
                    assignin('base',data.inputname,data.cst);
                end
                
                function editExponentAddCallback(hObj, event)
                    tag= str2num(get(hObj,'Tag')); %Nummer des VOIs
                    data.cst{tag,6}(size(data.cst{tag,6},2)).exponent = str2num(get(hObj, 'String'));
                    guidata(gcf, data);
                    % overwrite cst struct directly in workspace
                    assignin('base',data.inputname,data.cst);
                end
    
            end
                    
            function editPenaltyAddCallback(hObj, event)
                tag= str2num(get(hObj,'Tag')); %Nummer des VOIs 
                data.cst{tag,6}(size(data.cst{tag,6},2)).parameter(1) = str2num(get(hObj, 'String'));
                guidata(gcf, data);
                % overwrite cst struct directly in workspace
                assignin('base',data.inputname,data.cst);
            end
            
        end
       
      l=l+1;  
      
    end

    function pushbuttonDeleteCallback(hObj, event)
        tag = str2num(get(hObj,'Tag'));%Nummer des VOI, Nummer der Bedingung, Nummer der Zeile in Figure
        s = size(data.cst{tag(1),6},2);
        
        % delete in in data.cst
        data.cst{tag(1),6}(tag(2)) = [];
        guidata(gcf, data);

        % overwrite cst struct directly in workspace
        assignin('base',data.inputname,data.cst);
        
        if isempty(data.cst{tag(1),6})
            data.cst{tag(1),3} = 'IGNORED';
            guidata(gcf, data);

            % overwrite cst struct directly in workspace
            assignin('base',data.inputname,data.cst);
            
            dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.03 0.89-1/s*tag(3)*0.03 0.25 0.03]);
        
        end
        
        %löschen der Anzeige
        dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.27 0.89-tag(3)*0.03 0.075 0.03]);
        
        dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.39 0.9-tag(3)*0.03 0.075 0.02]);
        
        dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.51 0.9-tag(3)*0.03 0.075 0.02]);
        
        dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.63 0.9-tag(3)*0.03 0.075 0.02]);
        
        dataSetText = uicontrol('Style', 'text','String', {}, 'units', 'normalized',...
            'BackgroundColor', [0.8 0.8 0.8], 'Position', [0.75 0.9-tag(3)*0.03 0.075 0.02]);
        
    end

end