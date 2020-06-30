classdef matRad_Widget <  handle 
   
    properties (GetAccess = public , SetAccess = protected)
        widgetHandle    %Holds parent widget handle
        handles = struct([]);   %Holds all handles in parent handle
    end
    
    events
        workspaceChanged
    end
    
    methods
        %CONSTRUCTOR
        function this = matRad_Widget(handleParent)
            this.widgetHandle = handleParent;           
            this.createLayout();           
            this.initialize();
            
            [env, ~] = matRad_getEnvironment();
            
            % only enable in matlab
            %strcmp(env,'MATLAB') && 
            if strcmp(env,'MATLAB') && strcmp(get(handleParent,'type'),'figure')
                set(this.widgetHandle,'ButtonDownFcn',@(src,hEvent) update(this));   
                set(this.widgetHandle,'KeyPressFcn',@(src,hEvent) update(this));   
            end
        end

        function set.handles(obj,handles)
            obj.handles = handles;
        end
        
        %INITIALIZE FUNCTION
        function this = initialize(this)            
        end
        
        function this = update(this)
        end
        
        %{
        function this = disable(this)            
        end
        
        function this = enable(this)
        end
        %}
        
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
                Message = {Message,ME.message};%{Message,ME.getReport(meType,'hyperlinks','off')};
            end
            
            if isfield(handles,'ErrorDlg')
                if ishandle(handles.ErrorDlg)
                    close(handles.ErrorDlg);
                end
            end
            matRad_cfg.dispError(Message);
            %handles.ErrorDlg = errordlg(Message);
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
                Message = {Message,ME.message};%{Message,ME.getReport(meType,'hyperlinks','off')};
            end
            matRad_cfg.dispWarning(Message);
%
%             if isfield(handles,'ErrorDlg')
%                 if ishandle(handles.ErrorDlg)
%                     close(handles.ErrorDlg);
%                 end
%             end
%             handles.ErrorDlg = errordlg(Message);
%
             this.handles = handles;         

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
    end
    
    %methods (Abstract)
    %    this = initialize(this);
    %end

    methods (Access = protected)
        
        %CREATE LAYOUT FUNCTION
        function this = createLayout(this, handleParent)
        end  
        
        %Helper function to create handles structure from tags (similar to
        %guihandles)
        function this = createHandles(this)
            %Iterate through all objects
            all_h = findall(this.widgetHandle);
            
            for i = 1:numel(all_h)
                this_h = all_h(i);
                if isprop(this_h, 'Tag')
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
    end
end

