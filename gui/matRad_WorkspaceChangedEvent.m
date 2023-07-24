classdef  matRad_WorkspaceChangedEvent < handle
    %matRad_WorkspaceChangedEvent 
    %   Base class to store event data (no subclass yet for compatability)
    
    properties
        changedVariables %Cell array with changed variable names
    end
    
    methods
        function obj = matRad_WorkspaceChangedEvent(varargin)
            %matRad_WorkspaceChangedEvent Construct the event
            %   varargin is the cell array of arguments (strings)
            %   containing the variable names
            
            %obj = obj@event.EventData();
            
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

