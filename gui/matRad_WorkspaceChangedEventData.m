classdef  matRad_WorkspaceChangedEventData < matRad_WorkspaceChangedEvent & event.EventData
    %matRad_WorkspaceChangedEventData 
    %   EventData subclass to store changed variables when widget changes
    %   the workspace

    methods
        function obj = matRad_WorkspaceChangedEventData(varargin)
            %matRad_WorkspaceChangedEvent Construct the event
            %   varargin is the cell array of arguments (strings)
            %   containing the variable names
            
            obj = obj@event.EventData();
            obj = obj@matRad_WorkspaceChangedEvent(varargin{:});
            
            
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

