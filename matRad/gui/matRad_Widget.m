classdef matRad_Widget <  handle 

    % matRad_Widget Main Class for GUI widget generation 
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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    properties (GetAccess = public , SetAccess = protected)
        widgetHandle            %Holds parent widget handle
        handles = struct([]);   %Holds all handles in parent handle
        isInUifigure;
    end
    
    properties
        updateLock = false;     %Property to lock updating of the widget
    end
        
    events
        %If the widget changes the workspace, this event should be emitted 
        %with notify(...). Other widgets can add listeners to update when 
        %this event is emitted
        workspaceChanged 
    end
    
    methods
        %CONSTRUCTOR
        function this = matRad_Widget(handleParent)
            this.widgetHandle = handleParent;
            
            %Create the layout
            try 
                this.createLayout();           
            catch ME
                showError(this,'Widget could not be created!',ME);
            end

            %Initialize the widget
            try
                this.initialize();
            catch ME
                showWarning(this,'Widget could not be initialized!',ME);
            end
            
            matRad_cfg = MatRad_Config.instance();
            
            % only enable in matlab
            %strcmp(env,'MATLAB') && 
            if matRad_cfg.isMatlab && strcmp(get(handleParent,'type'),'figure')
                set(this.widgetHandle,'ButtonDownFcn',@(src,hEvent) update(this));   
                set(this.widgetHandle,'KeyPressFcn',@(src,hEvent) update(this));   
            end
        end

        function set.handles(obj,handles)
            obj.handles = handles;
        end
        
        %INITIALIZE FUNCTION
        function initialize(this)
            this.update();
        end
        
        function changedWorkspace(this,varargin)
            matRad_cfg = MatRad_Config.instance();
            % handle environment
            switch matRad_cfg.env
                case 'MATLAB'
                    %the PlanWidget only changes the pln
                    evt = matRad_WorkspaceChangedEventData(varargin{:});
                    notify(this, 'workspaceChanged',evt);
                case 'OCTAVE'
                    evt = matRad_WorkspaceChangedEvent(varargin{:});
                    matRad_notifyOctave(this, 'workspaceChanged',evt);
            end 
        end
        
        function handles = showError(this,Message,ME)
            matRad_cfg = MatRad_Config.instance();
            handles = this.handles;
            if nargin == 3
                %Add exception message
                if isfield(handles,'devMode') && handles.devMode
                    meType = 'extended';
                else
                    meType = 'basic';
                end
                Message = [Message,ME.message];%{Message,ME.getReport(meType,'hyperlinks','off')};
            end
            
            if isfield(handles,'ErrorDlg')
                if ishandle(handles.ErrorDlg)
                    close(handles.ErrorDlg);
                end
            end
            if ~matRad_cfg.disableGUI
                h = errordlg(Message);
                matRad_applyThemeToDlg(h);
            end
            matRad_cfg.dispError(Message);
            this.handles = handles;
        end
        
        function showWarning(this,Message,ME)
            matRad_cfg = MatRad_Config.instance();
            
            handles = this.handles;
            if nargin == 3
                %Add exception message
                if isfield(handles,'devMode') && handles.devMode
                    meType = 'extended';
                else
                    meType = 'basic';
                end
                Message = [Message,ME.message];
                % Future error hyperlinks {Message,ME.getReport(meType,'hyperlinks','off')};
            end
            if ~matRad_cfg.disableGUI
                h = warndlg(Message);
                matRad_applyThemeToDlg(h);
            end
            matRad_cfg.dispWarning(Message);
            this.handles = handles;         
        end

        function showMessage(this,message,varargin)
            matRad_cfg = MatRad_Config.instance();
            if ~matRad_cfg.disableGUI
                h = msgbox(message,varargin{:});
                matRad_applyThemeToDlg(h);
            end
        end
        
        %function notifyUpdate(this,workSpaceVariables)
        %   notify(this,'workspaceChanged',workSpaceVariables); 
        %end
        
        function enableWindowCallback(this,enable)
            
            if strcmp(get(this.widgetHandle,'type'),'figure') && enable
                set(this.widgetHandle,'ButtonDownFcn',@(src,hEvent) update(this));   
                set(this.widgetHandle,'KeyPressFcn',@(src,hEvent) update(this));
            else
                set(this.widgetHandle,'ButtonDownFcn','');   
                set(this.widgetHandle,'KeyPressFcn','');
            end
        end

        function isInUi = get.isInUifigure(this)
            hFig = ancestor(this.widgetHandle,'Figure');
            isInUi = matRad_isUifigure(hFig);
        end

        function delete(this)
            if ishandle(this.widgetHandle)
                this.destroyLayout();
            end        
        end

        function close(this)
            if ishandle(this.widgetHandle)
                this.destroyLayout();
            end        
        end
    end

    methods (Sealed)
        
        %Update function to be called when the widget should be updated
        %(i.e. by listener to the workspaceChanged event or in
        %initialization)
        %Implement the protected doUpdate function in a subclass to define
        %the update process
        function update(this,evt)
            try
                if nargin == 2
                    this.doUpdate(evt);
                else
                    this.doUpdate();
                end
            catch ME
                %Error Handling
                this.showWarning('Update failed',ME); %We only show a warning to not halt execution.
            end
        end
    end
    
    methods (Access = protected)
        
        function this = doUpdate(this,~)

        end

        %CREATE LAYOUT FUNCTION
        function this = createLayout(this, handleParent)
        end  

        function destroyLayout(this)
            fNames = fieldnames(this.handles);
            for i = 1:numel(fNames)
                h = this.handles.(fNames{i});
                if isscalar(h) && isgraphics(h)
                    delete(h);
                end
            end
        end
        
        %Helper function to create handles structure from tags (similar to
        %guihandles)
        function this = createHandles(this)
            %Iterate through all objects
            all_h = findall(this.widgetHandle);
            
            for i = 1:numel(all_h)
                this_h = all_h(i);
                if matRad_ispropCompat(this_h, 'Tag')
                    tag = get(this_h, 'Tag');
                    if ~isempty(tag) && isvarname(tag) % can it be used as a fieldname?
                        
                        % if a field of this name already exists, get its contents
                        if isfield(this.handles, tag)
                            prev_h = this.handles.(tag);
                        else
                            prev_h = [];
                        end
                        
                        % append our handle to whatever was there before. If nothing
                        % was there before.
                        if isappdata(this_h, 'Control')
                            % if this uicontrol is a proxy for external controls, replace it with
                            % that control
                            control = getappdata(this_h, 'Control');
                            this.handles(1).(tag) = [prev_h control.Instance];
                        else
                            this.handles(1).(tag) = [prev_h this_h];
                        end
                        
                    end % if legal tag
                end
            end % loop
        end
        
        %Helper function to check if update is necessary when an event is
        %passed
        function [doUpdate,whichVars] = checkUpdateNecessary(this,updateVars,evt)
            if isa(evt,'matRad_WorkspaceChangedEvent')
                changed = evt.changedVariables;
            else 
                changed = {};
            end
            
            
            %If changed is empty, we need an update, if not, we crosscheck
            doUpdate = true;
            whichVars = [];
                        
            if ~isempty(changed)
               cmp = cellfun(@(str) any(strcmp(str,changed)),updateVars);
               doUpdate = any(cmp);
               if nargout == 2
                   whichVars = find(cmp); %indices of variables
               end
            end            
        end
        
        %Helper function to position GUI elements on the Main Widget 
        function pos = computeGridPos(this,gridPos,buttonGridSize,buttonRelSize)
            if nargin < 4
                buttonRelSize = [0.9 0.75];
            end
        
            gridElementSize = 1./buttonGridSize;
            buttonRelSize = buttonRelSize ./ buttonGridSize;
            
            buttonGridPos = (gridPos-1) ./ buttonGridSize;
            buttonGridPos(2) = 1 - gridElementSize(2) - buttonGridPos(2);
            offset = (gridElementSize - buttonRelSize)./2;
            
            pos = [buttonGridPos+offset buttonRelSize];
        end
    end
    
end

