classdef  matRad_WorkspaceChangedEventData < matRad_WorkspaceChangedEvent & event.EventData
%  Class matRad_WorkspaceChangedEventData 
% EventData subclass to store changed variables when widget changes
% the workspace
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

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

