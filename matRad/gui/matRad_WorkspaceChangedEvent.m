classdef  matRad_WorkspaceChangedEvent < handle
% matRad_WorkspaceChangedEvent class 
% Base class to store event data (no subclass yet for compatability)
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

