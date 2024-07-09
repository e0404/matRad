function retval = matRad_addListenerOctave (hSource,eventName,callback)
% matRad_addListenerOctave is a function that creates and accumulates new
% Listeners in eventMap
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
    
    persistent warningprinted;
    if isempty(warningprinted)
        hSource.showWarning('You are using an experimental implementation of ''addlistener'' concept for Octave');
        warningprinted=2;
    end
    
    newListener = struct('src',hSource,'event',eventName,'callback',callback);
    if isempty(eventMap)
        eventMap = newListener;
    else
        eventMap(end+1) = newListener;
    end
    retval=newListener;
    
end
