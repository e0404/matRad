classdef matRad_Widget <  handle 
   
    properties (GetAccess = public , SetAccess = protected)
        widgetHandle    %Holds parent widget handle
        handles = [];   %Holds all handles in parent handle
    end
    
    methods
        %CONSTRUCTOR
        function this = matRad_Widget(handleParent)
            this.widgetHandle = handleParent;           
            this.createLayout();           
            this.initialize();
        end

        function set.handles(obj,handles)
            obj.handles = handles;
        end
        
        %INITIALIZE FUNCTION
        function this = initialize(this)
        end
    end

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
                            this.handles.(tag) = [prev_h control.Instance];
                        else
                            this.handles.(tag) = [prev_h this_h];
                        end
                        
                    end % if legal tag
                end
            end % loop
        end
    end
end

