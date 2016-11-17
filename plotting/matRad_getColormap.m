function output = matRad_getColormap(name,size)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function wrapper for getting a colormap. We use this wrapper to
% manually handle the supported colormaps enabling the definition of custom
% colormaps.
%
% call
%   cMap = matRad_getColormap(name,size)
%   list = matRad_getColormap()
% input
%   name        name of the colorbar
%   size        optional argument for the size / resolution of the colorbar
%   
%   if no argument is passed, a list (cell array) of the names of all 
%   supported colormaps will be returned
%
% output
%   This is either the requested colormap, or a list of all available
%   colormaps (see above)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    size = 64; %Default number of colors
end

%If no argument is given, the functinon returns all available colormaps
if nargin == 0
    %String cell array of supported colormaps
    cMaps{1}     = 'bone';
    cMaps{end+1} = 'gray';
    cMaps{end+1} = 'jet';
    cMaps{end+1} = 'parula';
    cMaps{end+1} = 'cool';
    
    % check for additional colormaps
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    listing = dir([folder filesep 'colormaps']);
    for i = 1:numel(listing)
        if listing(i).bytes > 0 && ~listing(i).isdir && strcmp(listing(i).name(end),'m')
            colormapFile = listing(i).name;
            cMaps{end+1} = colormapFile(1:end-2);
        end
    end
    output = cMaps;
end

if nargin > 0
    try
        sFuntion = [name '(' num2str(size,'%d') ')'];
        output = evalin('caller',sFuntion);
    catch
        warning(['Colormap "' name '" not supported by matRad']);
       output = jet(size);  
    end
end
