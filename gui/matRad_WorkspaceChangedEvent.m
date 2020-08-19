classdef  matRad_WorkspaceChangedEvent < event.EventData
    %matRad_WorkspaceChangedEvent 
    %   EventData subclass to store changed variables when widget changes
    %   the workspace
    
    properties
        changedVariables %Cell array with changed variable names
    end
    
    methods
        function obj = matRad_WorkspaceChangedEvent(varargin)
            %matRad_WorkspaceChangedEvent Construct the event
            %   varargin is the cell array of arguments (strings)
            %   containing the variable names
            
            obj = obj@event.EventData();
            
            obj.changedVariables = varargin;
            
            
            %Debug:
            %{
            if isempty(varargin)
                changed = 'all';
            else
                changed = strjoin(varargin,'|');
            end
            fprintf('Changed variables: %s\n',changed);
            %}
        end

    end
end

