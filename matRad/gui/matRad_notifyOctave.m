function matRad_notifyOctave(hObject,eventName,evt)
% Experimental Function to notify on event for OCTAVE
%
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

    global eventMap;
    if isempty(eventMap)
        return;
    end
        
    persistent warningprinted;
    if isempty(warningprinted)
        hObject.showWarning('You are using an experimental implementation of ''notify'' concept for Octave');
        warningprinted=2;
    end
    
    objHandles = {eventMap(:).src};
    
    %clean Handles
    validH = cellfun(@isa,objHandles,repmat({'handle'},size(objHandles))); %cellfun(@ishandle,objHandles);
    objHandles = objHandles(validH);
    eventMap = eventMap(validH);
    
    %Find the object    
    objIx = cellfun(@(h) isequal(hObject,h), objHandles);
    
    objEvents = eventMap(objIx);
    
    if isempty(objEvents)
        return;
    end
    
    allNames = {objEvents(:).event};
    eventIx = find(strcmp(eventName,allNames)); 
    
    for runIx = 1:numel(eventIx)
        runEventIx = eventIx(runIx);
        runEvent = objEvents(runEventIx);
        if nargin < 3
            runEvent.callback(hObject);
        else
            runEvent.callback(hObject,evt);
    end 
end
